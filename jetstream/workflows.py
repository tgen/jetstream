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
import os
import re
import logging
import shutil
from datetime import datetime
import networkx as nx
from networkx.readwrite import json_graph
import jetstream
from threading import Lock
from collections import Counter, deque
from jetstream import utils

log = logging.getLogger(__name__)


class NoFallback:
    pass


def search_pattern(pat):
    return re.compile('^{}$'.format(pat))


def save(workflow, path):
    start = datetime.now()
    lock_path = path + '.lock'

    with open(lock_path, 'w') as fp:
        log.critical('Saving workflow...'.format(lock_path))
        fp.write(workflow.to_yaml())

    shutil.move(lock_path, path)

    elapsed = datetime.now() - start
    log.critical('Elapsed: {} Workflow saved to {}'.format(elapsed, path))


def coerce_str_to_list(value):
    if isinstance(value, str):
        return [value, ]
    else:
        return value


def coerce_str_to_bytes(value):
    if isinstance(value, str):
        return value.encode('utf8')
    else:
        return value


class Task(object):
    states_lookup = {
        0: 'new',
        1: 'pending',
        2: 'complete',
        3: 'failed'
    }

    valid_states = list(states_lookup.keys())
    status_lookup = {v: k for k, v in states_lookup.items()}
    valid_status = list(status_lookup.keys())

    def __init__(self, id, *, cmd=None, before=None, after=None, input=None,
                 output=None, stdin=None, stdout=None, stderr=None, cpus=0,
                 mem=0, walltime=0, status='new', returncode=None, start=None,
                 end=None, methods=None, description=None, help=None):

        self.id = str(id)
        self.cmd = cmd
        self.before = coerce_str_to_list(before)
        self.after = coerce_str_to_list(after)
        self.input = coerce_str_to_list(input)
        self.output = coerce_str_to_list(output)
        self.stdin = stdin
        self.stdout = stdout

        if self.stdout is None:
            # Remove whitespace from id to form default out path
            base_path = os.path.join('logs', re.sub(r'\s', '_', self.id))
            self.stdout = base_path + '.out'

        self.stderr = stderr or self.stdout
        self.cpus = int(cpus)
        self.mem = mem
        self._state = None
        self.status = status
        self.returncode = returncode
        self.start = start
        self.end = end
        self.walltime = walltime
        self.methods = methods
        self.description = description
        self.help = help
        self._workflow = None

    def __repr__(self):
        return '<Task({}): {} >'.format(self.status, self.id)

    def __hash__(self):
        return hash(self.id)

    def pretty(self):
        return utils.yaml_dumps(self.serialize())

    def serialize(self):
        data = {k: v for k, v in vars(self).items() if not k.startswith('_')}
        data.update(status=self.status)
        return data

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, value):
        if value not in Task.states_lookup:
            err = 'State must be one of {}'.format(Task.valid_states)
            raise ValueError(err) from None
        else:
            self._state = value

    @property
    def status(self):
        return Task.states_lookup[self._state]

    @status.setter
    def status(self, value):
        try:
            self._state = Task.status_lookup[value]
        except KeyError:
            err = 'Status must be one of {}'.format(Task.valid_status)
            raise KeyError(err) from None

    def is_new(self):
        if self.status == 'new':
            return True
        else:
            return False

    def is_pending(self):
        if self.status == 'pending':
            return True
        else:
            return False

    def is_complete(self):
        if self.status in ('complete', 'failed'):
            return True
        else:
            return False

    def is_ready(self):
        return self._workflow.is_ready(self.id)

    def reset(self):
        log.critical('{} reset!'.format(self.id))

        self._state = 0
        self.returncode = None
        self.start = None
        self.end = None

    def pending(self):
        log.critical('{} is pending!'.format(self.id))

        self._state = 1
        self.start = str(datetime.now())

    def complete(self, returncode=0):
        log.critical('{} is complete!'.format(self.id))

        self._state = 2
        self.returncode = returncode
        self.end = str(datetime.now())

    def fail(self, returncode=1):
        log.critical('{} failed!'.format(self.id))

        self._state = 3
        self.returncode = returncode
        self.end = str(datetime.now())

        for task in self._workflow.dependents(self.id):
            task.fail(returncode)

    def add_to_workflow(self, wf):
        self._workflow = wf


class Workflow(object):
    def __init__(self):
        self.graph = nx.DiGraph()
        self._lock = Lock()
        self._stack = list()

    def __enter__(self):
        """ Workflows can be edited in a transaction using the context manager
        statement "with". This allows multiple task additions to take place
        with only a single update to the workflow edges. """
        self._lock.acquire()
        self._stack = list()

    def __exit__(self, exc_type, exc_value, traceback):
        """ If there is an error during the transaction, all nodes will
        be rolled back on exit. """

        if exc_value is not None:
            for task_id in self._stack:
                del self.graph.nodes()[task_id]

        self.update()
        self._stack = list()
        self._lock.release()

    def __repr__(self):
        stats = Counter([t.status for t in self.tasks(objs=True)])
        return '<jetstream.Workflow {}>'.format(stats)

    def __len__(self):
        return len(self.graph)

    def __iter__(self):
        return WorkflowIterator(self)

    def add_task(self, task_id, **directives):
        if task_id in self.graph:
             raise ValueError('Duplicate task id: {}'.format(task_id))

        t = Task(task_id, **directives)
        t.add_to_workflow(self)

        self.graph.add_node(task_id, obj=t)

        if self.is_locked():
            self._stack.append(task_id)
        else:
            try:
                self.update()
            except Exception as e:
                self.graph.remove_node(task_id)
                raise e

        return t

    def update(self):
        current = list(self.graph.edges())
        self.graph.remove_edges_from(current)

        try:
            for task_id in self.graph.nodes():
                self._link_dependencies(task_id)
        except Exception as e:
            self.graph.remove_edges_from(list(self.graph.edges()))
            self.graph.add_edges_from(current)
            raise e from None

    def resume(self):
        """ Returns all pending nodes to an incomplete state. """
        for task in self.tasks(objs=True):
            if task.status == 'pending':
                task.reset()

    def reset(self):
        """ Returns all nodes to a new state. """
        for task in self.tasks(objs=True):
            task.reset()

    def retry(self):
        """ Resets all pending and failed tasks. """
        for task in self.tasks(objs=True):
            if task.status in ('pending', 'failed'):
                task.reset()

    def tasks(self, objs=False):
        if objs:
            return (t['obj'] for i, t in self.graph.nodes(data=True))
        else:
            return self.graph.nodes()

    def get_task(self, task_id):
        return self.graph.nodes[task_id]['obj']

    def dependencies(self, task_id):
        return (self.get_task(dep) for dep in self.graph.successors(task_id))

    def dependents(self, task_id):
        return (self.get_task(dep) for dep in self.graph.predecessors(task_id))

    def is_ready(self, task_id):
        """ Returns True if "task_id" is ready for execution. """
        task = self.get_task(task_id)

        if task.status != 'new':
            return False

        for dependency in self.dependencies(task_id):
            if dependency.status != 'complete':
                return False
        else:
            return True

    def find_by_id(self, pattern, fallback=NoFallback):
        log.debug('Find by id pattern: {}'.format(pattern))

        pat = search_pattern(pattern)
        fn = lambda task_id: pat.match(task_id)
        gen = self.graph.nodes()
        matches = set(filter(fn, gen))

        log.debug('Found matches: {}'.format(matches))

        if matches:
            return matches
        elif fallback is NoFallback:
            raise ValueError('No tasks match value: {}'.format(pattern))
        else:
            return fallback

    def find_by_cmd(self, pattern, fallback=NoFallback):
        log.debug('Find by cmd pattern: {}'.format(pattern))

        pat = search_pattern(pattern)
        matches = set()

        for task_id, data in self.graph.nodes(True):
            task = data['obj']

            if task.cmd is None:
                continue

            if pat.match(task.cmd):
                matches.add(task_id)

        if matches:
            return matches
        elif fallback is NoFallback:
            raise ValueError('No tasks match value: {}'.format(pattern))
        else:
            return fallback

    def find_by_output(self, pattern, fallback=NoFallback):
        log.debug('Find by output pattern: {}'.format(pattern))

        pat = search_pattern(pattern)
        matches = set()

        for task_id, data in self.graph.nodes(True):
            task = data['obj']

            if task.output is None:
                continue

            for value in task.output:
                if pat.match(value):
                    matches.add(task_id)

        if matches:
            return matches
        elif fallback is NoFallback:
            raise ValueError('No tasks match value: {}'.format(pattern))
        else:
            return fallback

    def is_locked(self):
        return self._lock.locked()

    def serialize(self):
        data = to_node_link_data(self)

        for node in data['nodes']:
            node['obj'] = node['obj'].serialize()

        return data

    def to_yaml(self):
        return utils.yaml_dumps(self.serialize())

    def dump_yaml(self, *args, **kwargs):
        return utils.yaml_dump(self.serialize(), *args, **kwargs)

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

    def _link_dependencies(self, task_id):
        log.debug('Linking dependencies for: {}'.format(task_id))
        task = self.get_task(task_id)

        if task.after:
            log.debug('Linking "after" dependencies')
            # "after" specifies edges that should run:
            #    task_id ---depends on---> target, ...

            for value in task.after:
                matches = self.find_by_id(value)

                if task.id in matches:
                    raise ValueError(
                        'Task "after" directives cannot match itself - '
                        'Task: {} Pattern: {}'.format(task_id, value)
                    )

                for tar_id in matches:
                    self._add_edge(task_id, tar_id)

        if task.before:
            log.debug('Linking "before" dependencies')

            # "before" specifies edges that should run:
            #    task_id <---depends on--- target, ...

            for value in task.before:
                matches = self.find_by_id(value)

                if task_id in matches:
                    raise ValueError(
                        'Task "before" directives cannot match itself - '
                        'Task: {} Pattern: {}'.format(task_id, value)
                    )

                for tar_id in matches:
                    self._add_edge(tar_id, task_id)

        if task.input:
            log.debug('Linking "input" dependencies')

            # "input" specifies edges that should run:
            #    task_id ---depends on---> target, ...
            # Where target includes an "output" value
            # matching the "input" value.

            for value in task.input:
                matches = self.find_by_output(value)

                if task_id in matches:
                    raise ValueError(
                        'Task "input" directives cannot match itself - '
                        'Task: {} Pattern: {}'.format(task_id, value)
                    )

                for tar_id in matches:
                    self._add_edge(task_id, tar_id)

    def compose(self, wf):
        """ Compose this workflow with another.

        This add nodes and edges from another workflow to this workflow
        graph. See networkx.compose(). G -> H where G is "new_wf" and H
        is this workflow. The result is that node status from the current
        workflow is preserved. If a workflow contains overlapping node ids
        this will essentially extend the workflow.

        ::

                 G (new_wf) -->    H (self)     =   self.graph
            ---------------------------------------------------
                (A)new                              (A)complete
                 |         --->   (A)complete   =    |
                (B)new                              (B)new


        :param new_wf: Another workflow to add to this workflow
        :return: None
        """
        old_graph = self.graph
        graph = nx.compose(wf.graph, self.graph)

        try:
            self.graph = graph
            self.update()
        except Exception as e:
            log.critical('Composition would result in a workflow with '
                         'errors in task dependency graph. See traceback '
                         'for more details.')
            self.graph = old_graph
            raise e from None

    def compose_all(self, *wfs):
        """ Compose this workflow with multiple other workflows.
        See Workflow.compose for more details. """
        old_graph = self.graph
        wfs = [wf.graph for wf in wfs] + [self.graph]
        graph = nx.compose_all(wfs)

        try:
            self.graph = graph
            self.update()
        except Exception as e:
            log.critical('Composition would result in a workflow with '
                         'errors in task dependency graph. See traceback '
                         'for more details.')
            self.graph = old_graph
            raise e from None

    def pretty(self):
        return utils.yaml_dumps(self.serialize())


class WorkflowIterator(object):
    def __init__(self, workflow):
        self.workflow = workflow
        self._queues = [deque(), deque()]
        self.complete = list()
        self.queue, self.stack = self._queues

        self.workflow.retry()

        for task in self.workflow.tasks(objs=True):
            if task.is_complete():
                self.complete.append(task)
            else:
                self.queue.append(task)

    def __repr__(self):
        return '<jetstream.Workflow {}>'.format(self.status())

    def status(self):
        return {
            'queue': len(self.queue),
            'stack': len(self.stack),
            'complete': len(self.complete)
        }

    def swap(self):
        log.debug('Swap queue and stack')
        self.queue, self.stack = self.stack, self.queue

    def __next__(self):
        there_are_pending_tasks = False
        checked_stack = False

        while 1:
            try:
                task = self.queue.popleft()
                log.debug('Checking {}'.format(task.id))

                if task.is_complete():
                    self.complete.append(task)

                elif task.is_pending():
                    self.stack.append(task)
                    there_are_pending_tasks = True

                elif task.is_ready():
                    self.stack.append(task)
                    task.pending()
                    return task

                else:
                    # task is not ready
                    self.stack.append(task)

            except IndexError:
                log.debug('Queue empty, checking stack...')

                if not self.stack:
                    if there_are_pending_tasks:
                        return None
                    else:
                        # Queue is empty, stack is empty, and no pending.
                        log.critical('All tasks complete!')
                        raise StopIteration from None

                elif self.stack and not checked_stack:
                    checked_stack = True
                    self.swap()

                else:
                    return None


def from_node_link_data(data):
    graph = json_graph.node_link_graph(data)
    wf = Workflow()

    for node_id, node_data in graph.nodes(data=True):
        task_id = node_data['obj'].pop('id')
        wf.add_task(task_id, **node_data['obj'])

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


def build_workflow(tasks):
    """ Given a sequence of tasks (dictionaries with properties described in
    the workflow specification), returns a workflow with nodes and edges
    already added """
    log.critical('Building workflow...')

    if isinstance(tasks, str):
        # If tasks are not parsed yet, do it automatically
        tasks = utils.yaml_loads(tasks)

    if not tasks:
        raise ValueError('No tasks were found in the data!')

    wf = Workflow()

    with wf:
        for task in tasks:
            wf.add_task(task_id=task.pop('id'), **task)

    log.critical('Workflow ready: {}'.format(wf))
    return wf
