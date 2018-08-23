"""Network graph model of computational workflows

The `Workflow` class models computational workflows as a directed-acyclic graph,
where nodes are tasks to complete, and edges represent dependencies between
those tasks. It includes methods for building workflows (add_task,
add_dependency) in addition to the methods required for executing a workflow
(__next__, fail, complete, reset, etc.).


Workflows are built by rendering a template with data. Templates are text
documents that describe a set of tasks to complete, and they can include
dynamic elements through templating syntax. Data can be saved in files located
in the config directory of project, or given as arguments to template (see
Project.render or Project.run).

..

    Template + Data --Render--> Workflow

Building a Workflow
--------------------

Templates are a set of tasks described in YAML format:

.. code-block:: yaml

   - id: align_fastqs
    cmd: bwa mem grch37.fa sampleA_R1_001.fastq.gz sampleA_R2_001.fastq.gz

   - id: index_bams
     cmd: samtools index sampleA.bam


Dependencies can be specified in a task with "before" or "after":

.. code-block:: yaml

  - id: align_fastqs
    cmd: bwa mem grch37.fa sampleA_R1_001.fastq.gz sampleA_R2_001.fastq.gz

  - id: index_bams
    cmd: samtools index sampleA.bam
    after: align_fastqs


or, equivalently:

.. code-block:: yaml

  - id: align_fastqs
    cmd: bwa mem grch37.fa sampleA_R1_001.fastq.gz sampleA_R2_001.fastq.gz
    before: index_bams

  - id: index_bams
    cmd: samtools index sampleA.bam


Finally, Jinja templating can be used to add dynamic elements.

.. code-block:: ansible

    {% for sample in project.config.samples %}

    - id: align_fastqs
      cmd: bwa mem ref.fa {{ sample.r1_fastq }} {{ sample.r2_fastq }}

    - id: index_bam_{{ sample.name }}
      cmd: samtools index {{ sample.name}}.bam

    {% endfor %}


After tasks are rendered with data, a Workflow is built from the tasks.
Workflows are built by adding a each task to the workflow, then adding
dependencies for every "before" or "after" directive found in the tasks.
Internally, this creates a directed-acyclic graph where nodes represent tasks
to complete, and edges represent dependencies between those tasks.


Upfront workflow rendering
---------------------------

Rendering a template is a dynamic procedure that uses input data and template
directives to generate a set of tasks. But, after those tasks are used to build
a workflow, the resulting workflow is a final, complete, description of the
commands required, and the order in which they should be executed.

Workflows do not change in response to events that occur during runtime. If a
task exists in a workflow, the runner will always launch it. Unlike other
worklow engines, there is no flow control (conditionals, loops, etc.) contained
in the tasks themselves. Flexibility is enabled by templates. The only
exception to this is that a task will automatically fail if any of its
dependencies fails.

Cases where flexibility during runtime may be necessary:

Some input data needs to be split into n chunks where n is determined
by a command during runtime. Each chunk then needs to be treated as an
individual task in the workflow. Note that this is not a problem if n can be
determined prior to runtime, or if the command can handle the chunking
internally.

"""
import re
import logging
import shutil
from datetime import datetime
import networkx as nx
from networkx.readwrite import json_graph
import jetstream
from threading import Lock
from collections import Counter
from jetstream import utils, log
from jetstream.tasks import Task

log = logging.getLogger(__name__)


class Workflow(object):
    def __init__(self):
        self.graph = nx.DiGraph()
        self._lock = Lock()
        self._stack = list()

    def __contains__(self, item):
        if isinstance(item, str):
            return item in self.graph
        else:
            return item.tid in self.graph

    def __enter__(self):
        """Workflows can be edited in a transaction using the context manager
        statement "with". This allows multiple task additions to take place
        with only a single update to the workflow edges. """
        self._lock.acquire()
        self._stack = list()

    def __exit__(self, exc_type, exc_value, traceback):
        """If there is an error during the transaction, all nodes will
        be rolled back on exit. """
        if exc_value is not None:
            for task_id in self._stack:
                self.graph.remove_node(task_id)

        self.update()
        self._stack = list()
        self._lock.release()

    def __iter__(self):
        return WorkflowIterator(self)

    def __len__(self):
        return len(self.graph)

    def __repr__(self):
        stats = Counter([t.status for t in self.tasks(objs=True)])
        return '<jetstream.Workflow {}>'.format(stats)

    def _add_edge(self, from_node, to_node):
        """ Edges represent dependencies between tasks. Edges run FROM one node
        TO another node that it depends upon. Nodes can have multiple edges,
        but not multiple instances of the same edge (multigraph).

            Child ----- Depends Upon -----> Parent
         (from_node)                       (to_node)

        This means that the out-degree of a node represents the number of
        dependencies it has. A node with zero out-edges is a "root" node, or a
        task with no dependencies.

        The "add_dependency" method is provided for adding dependencies to a
        workflow, and should be preferred over adding edges directly to the
        workflow graph. """
        log.debug('Adding edge: {} -> {}'.format(from_node, to_node))

        self.graph.add_edge(from_node, to_node)

        if not nx.is_directed_acyclic_graph(self.graph):
            self.graph.remove_edge(from_node, to_node)
            raise jetstream.NotDagError

    def _add_node(self, task):
        """Add a node to the graph.
        Nodes are expected to be an instance of Jetstream.Task. If the workflow
        is not locked, via "with" statement, this will trigger Workflow.update.
        """
        if not isinstance(task, Task):
            raise ValueError('task must be instance of {}'.format(Task))

        if task.tid in self.graph:
            raise ValueError('Duplicate task ID: {}'.format(task.tid))

        log.info('Adding task: {}'.format(task))

        task.workflow = self
        self.graph.add_node(task.tid, obj=task)

        if self.is_locked():
            self._stack.append(task.tid)
        else:
            try:
                self.update()
            except Exception as e:
                self.graph.remove_node(task.tid)
                raise e

        return task

    def _flow_after(self, task):
        """Generate edges for "after" directives of a task
        "after" specifies edges that run:
            task ---depends on---> target, ...
        """
        after = task.get('after')
        log.debug('"after" directive: {}'.format(after))

        if after:
            after = utils.coerce_sequence(after)
            log.debug('"after" directive after coercion: {}'.format(after))

            for value in after:
                matches = self.find(value)

                if task.tid in matches:
                    raise ValueError(
                        'Task "after" directives cannot match itself - '
                        'Task: {} Pattern: {}'.format(task.tid, value)
                    )

                for tar_id in matches:
                    self._add_edge(task.tid, tar_id)

    def _flow_before(self, task):
        """Generate edges for "before" directives of a task
        "before" specifies edges that should run:
            task <---depends on--- target, ...
        """
        before = task.get('before')
        log.debug('"before" directive: {}'.format(before))

        if before:
            before = utils.coerce_sequence(before)
            log.debug('"before" directive after coercion: {}'.format(before))

            for value in before:
                matches = self.find(value)

                if task.tid in matches:
                    raise ValueError(
                        'Task "before" directives cannot match itself - '
                        'Task: {} Pattern: {}'.format(task.tid, value)
                    )

                for tar_id in matches:
                    self._add_edge(tar_id, task.tid)

    def _flow_input(self, task):
        """Generate edges for "input" directives of a task
        "input" specifies edges that should run:
            task ---depends on---> target, ...
        Where target includes an "output" value matching the "input" value."""
        input = task.get('input')
        log.debug('"input" directive: {}'.format(input))

        if input:
            input = utils.coerce_sequence(input)
            log.debug('"input" directive after coercion: {}'.format(input))

            for value in input:
                matches = self.find_by_output(value)

                if task.tid in matches:
                    raise ValueError(
                        'Task "input" directives cannot match itself - '
                        'Task: {} Pattern: {}'.format(task.tid, value)
                    )

                for tar_id in matches:
                    self._add_edge(task.tid, tar_id)

    def add_task(self, task):
        return self._add_node(task)

    def compose(self, workflow):
        """Compose this workflow with another workflow.
        This adds all tasks from another workflow to this workflow.

        Note: This will not roll back changes if an error occurs during
        the merge.

        ::

                G (wf)    --->    H (self)     =   self.graph
           ---------------------------------------------------
                                    (A)new         (A)complete
             (A)complete  --->       |        =     |
                                    (B)new         (B)new

        :param workflow: Another workflow to add to this workflow
        :return: None
        """
        log.info('Composing {} with {}'.format(self, workflow))

        with self:
            for task in workflow.tasks(objs=True):
                if not task in self:
                    self.add_task(task)
                else:
                    log.info('{} already in workflow'.format(task))

    def dependencies(self, task):
        """Returns a generator that yields the dependencies of a given task"""
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.tid

        return (self.get_task(tid) for tid in self.graph.successors(task_id))

    def dependents(self, task):
        """Returns a generator that yields the dependents of a given task"""
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.tid

        return (self.get_task(tid) for tid in self.graph.predecessors(task_id))

    def find(self, pattern, fallback=utils.sentinel):
        """Find searches for tasks by "name" with a regex pattern."""
        log.debug('Find: {}'.format(pattern))

        pat = search_pattern(pattern)
        matches = set()

        for task_id, data in self.graph.nodes(True):
            task = data['obj']
            name = task.get('name')

            if name is None:
                continue

            if pat.match(name):
                matches.add(task_id)

        if matches:
            return matches
        elif fallback is utils.sentinel:
            raise ValueError('No task names match value: {}'.format(pattern))
        else:
            return fallback

    # TODO I've accumulated several "find" methods that can probably be
    # generalized into a find_by(key, value) method.

    def find_by_id(self, pattern, fallback=utils.sentinel):
        log.debug('Find by id pattern: {}'.format(pattern))

        pat = search_pattern(pattern)
        fn = lambda task_id: pat.match(task_id)
        gen = self.graph.nodes()
        matches = set(filter(fn, gen))

        log.debug('Found matches: {}'.format(matches))

        if matches:
            return matches
        elif fallback is utils.sentinel:
            raise ValueError('No task ids match value: {}'.format(pattern))
        else:
            return fallback

    def find_by_output(self, pattern, fallback=utils.sentinel):
        log.debug('Find by output: {}'.format(pattern))

        pat = search_pattern(pattern)
        matches = set()

        for task_id, data in self.graph.nodes(True):
            task = data['obj']
            output = task.get('output')

            if output is None:
                continue

            output = utils.coerce_sequence(output)

            for value in output:
                if pat.match(value):
                    matches.add(task_id)
                    break

        if matches:
            return matches
        elif fallback is utils.sentinel:
            raise ValueError('No task outputs match value: {}'.format(pattern))
        else:
            return fallback

    def get_task(self, task_id):
        return self.graph.nodes[task_id]['obj']

    def is_locked(self):
        """Returns True if this workflow is locked.
        A locked workflow will not automatically trigger edge updates when
        adding new tasks. This allows tasks to be added "out-of-order"
        regarding their dependency chain. """
        return self._lock.locked()

    def is_ready(self, task):
        """Returns True if task is ready for execution."""
        if isinstance(task, str):
            task = self.get_task(task)

        if task.status != 'new':
            return False

        for dependency in self.dependencies(task):
            if not dependency.is_done():
                return False
        else:
            return True

    def list_tasks(self):
        return list(self.tasks(objs=True))

    def new_task(self, *args, **kwargs):
        """Shortcut to create a new Task object and add to this workflow."""
        task = Task(*args, **kwargs)
        return self.add_task(task)

    def remove_task(self, task_id):
        return self.graph.remove_node(task_id)

    def reset(self):
        """Resets all tasks state."""
        log.critical('Resetting state for all tasks')
        for task in self.tasks(objs=True):
            task.reset()

    def resume(self):
        """Resets all "pending" tasks state."""
        for task in self.tasks(objs=True):
            if task.status == 'pending':
                task.reset()

    def retry(self):
        """Resets all "pending" and "failed" tasks state."""
        for task in self.tasks(objs=True):
            if task.status in ('pending', 'failed'):
                task.reset()

    def serialize(self):
        """Convert the workflow to a node-link formatted object that can
        be easily dumped to JSON/YAML """
        data = to_node_link_data(self)

        for node in data['nodes']:
            node['obj'] = node['obj'].serialize()

        return data

    def tasks(self, objs=True):
        """Access to the tasks in this workflow.
        If objs is False, only the tids will be returned."""
        if objs:
            return (t['obj'] for i, t in self.graph.nodes(data=True))
        else:
            return self.graph.nodes()

    def to_json(self, indent=None, sort_keys=True):
        """Returns JSON string representation of this workflow"""
        return utils.json_dumps(
            self.serialize(), indent=indent, sort_keys=sort_keys
        )

    def to_yaml(self):
        """Returns YAML string representation of this workflow"""
        return utils.yaml_dumps(self.serialize())

    def update(self):
        """Recalculate the edges for this workflow"""
        for task in self.tasks(objs=True):
            log.debug('Linking dependencies for: {}'.format(task))
            self._flow_after(task)
            self._flow_before(task)
            self._flow_input(task)


class WorkflowIterator(object):
    def __init__(self, workflow):
        self.workflow = workflow
        self.total = len(workflow)
        self.tasks = workflow.list_tasks()
        self.pending = list()

    def __repr__(self):
        msg = '{}/{} remaining'.format(
            len(self.tasks) + len(self.pending), self.total)
        return '<jetstream.WorkflowIterator: {}>'.format(msg)
    
    def __next__(self):
        log.verbose('Request for next task')
        
        self.pending = [t for t in self.pending if not t.is_done()]
        log.verbose('Pending: {}'.format(self.pending))
        
        if not self.tasks and not self.pending:
            raise StopIteration
        
        for i in reversed(range(len(self.tasks))):
            task = self.tasks[i]
            log.verbose('Considering: {}'.format(task))
            
            if task.is_done():
                self.tasks.pop(i)
                log.info(self.workflow)
            elif task.is_ready():
                self.tasks.pop(i)
                self.pending.append(task)
                task.start()
                return task
        else:
            return None


def search_pattern(pat):
    return re.compile('^{}$'.format(pat))


def save_workflow(workflow, path):
    start = datetime.now()
    lock_path = path + '.lock'

    with open(lock_path, 'w') as fp:
        log.info('Saving workflow...'.format(lock_path))
        fp.write(workflow.to_yaml())

    shutil.move(lock_path, path)

    elapsed = datetime.now() - start
    log.info('Workflow saved (after {}) to {}'.format(elapsed, path))


def from_node_link_data(data):
    graph = json_graph.node_link_graph(data)
    wf = Workflow()

    for node_id, node_data in graph.nodes(data=True):
        wf.new_task(node_data['obj'])

    return wf


def to_node_link_data(wf):
    return json_graph.node_link_data(wf.graph)


def to_cytoscape_json_data(wf):
    """ Export a workflow as a cytoscape JSON file

    Cytoscape is good for vizualizing network graphs. It complains about node
    data that are not strings, so all node data are converted to strings on
    export. This causes Cytoscape json files to be a one-way export, they
    cannot be loaded back as workflow. """
    data = nx.cytoscape_data(wf.graph)

    for n in data['elements']['nodes']:
        for k, v in n['data'].items():
            n['data'][k] = str(v)

    return data


def load_workflow(path):
    """ Load a workflow from a file. """
    graph = utils.yaml_load(path)
    return from_node_link_data(graph)


def build_workflow_from_string(tasks):
    parsed = utils.yaml_loads(tasks)
    return build_workflow(parsed)


def build_workflow(tasks):
    """ Given a sequence of tasks (dictionaries with properties described in
    the workflow specification), returns a workflow with nodes and edges
    already added """
    log.info('Building workflow...')

    if isinstance(tasks, str):
        # If tasks are not parsed yet, do it automatically
        tasks = utils.yaml_loads(tasks)

    if not tasks:
        raise ValueError('No tasks were found in the data!')

    wf = Workflow()

    with wf:
        for task_mapping in tasks:
            wf.add_task(Task(data=task_mapping))

    log.info('Workflow ready: {}'.format(wf))
    return wf
