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
import time
import logging
import shutil
from datetime import datetime
import networkx as nx
from networkx.readwrite import json_graph
import jetstream
from collections import Counter
from jetstream import utils

log = logging.getLogger(__name__)


def save(workflow, path):
    start = datetime.now()
    lock_path = path + '.lock'
    data = to_node_link_data(workflow)

    with open(lock_path, 'w') as fp:
        log.critical('Saving workflow...'.format(lock_path))
        utils.yaml.dump(data, fp)

    shutil.move(lock_path, path)

    elapsed = datetime.now() - start
    log.critical('Elapsed: {} Workflow saved to {}'.format(elapsed, path))


def autosave(f):
    """ Decorator for workflow methods that should cause the workflow
    be saved. """
    err = 'Autosave requires Workflow.path or Workflow.project to be set'

    def decorator(workflow, *args, **kwargs):
        if workflow.autosave:
            if workflow.path is None and workflow.project is None:
                raise ValueError(err)

            now = time.time()

            if (now - workflow._last_save) > workflow.save_interval:
                res = f(workflow, *args, **kwargs)
                save(workflow, workflow.path or workflow.project.workflow_path)
                workflow._last_save = time.time()
                return res
            else:
                return f(workflow, *args, **kwargs)

        else:
            return f(workflow, *args, **kwargs)

    return decorator


class Workflow:
    """ Workflows are a network graph representing a series of commands that
    need to be executed on a project. Each node of the graph represents a task
    to be completed, and edges represent dependencies between the tasks.

    Workflows are iterables that implement a ``__next__`` method. Calling
    ``next()`` on a workflow will return the next task, or None if there are
    none currently available. StopIteration will be raised when all tasks are
    complete.

    Task status should be updated via Workflow.fail, or Workflow.complete. This
    will ensure any dependencies are updated accordingly. Workflow.update can
    be used update other task properties.

    Workflows can be loaded from existing graphs with the data argument or
    from_node_link_data() method. """
    def __init__(self, project=None, graph=None, path=None, autosave=False,
                 save_interval=300):
        self.project = project
        self.graph = graph or nx.DiGraph(_backup=dict())
        self.path = path
        self.autosave = autosave
        self.save_interval = save_interval
        self._new_tasks = None
        self._pending_tasks = None
        self._complete_tasks = None
        self._last_save = time.time()

        if not nx.is_directed_acyclic_graph(self.graph):
            raise jetstream.NotDagError

    def __repr__(self):
        return '<jetstream.Workflow {}>'.format(Counter(self.status()))

    def __len__(self):
        return len(self.graph)

    def __iter__(self):
        return WorkflowIterator(self)

    def _root_nodes(self):
        """ Returns the set of root nodes in the graph. """
        res = set()
        for node in self.graph.nodes():

            # Note: code inspectors warn about calling graph.out_degree
            # because (I think) it's an overloaded method. But this line seems
            # to work as intended. The goal is to find nodes in the graph with
            # an out_degree of 0.

            if self.graph.out_degree(node) == 0:
                res.add(node)

        return res

    def _backup_node(self, node_id, data):
        task_data = tuple(sorted(data.copy().items()))
        self.graph.graph['_backup'][node_id] = task_data

    def _restore_node(self, node_id):
        try:
            data = dict(self.graph.graph['_backup'][node_id])
            self.graph.nodes[node_id].clear()
            self.graph.nodes[node_id].update(data)
        except KeyError:
            log.critical('Unable to restore node from backup.')
            self.graph.nodes[node_id].update(status='new')

    def _add_node(self, node_id, data):
        """ Adding a node requires unique node_id. A RuntimeError will be
        raised if the node_id already exists in the graph. """
        log.debug('Adding node: {}'.format(node_id))

        if self.get_task(node_id):
            raise ValueError('Duplicate node id: {}'.format(node_id))

        # All nodes get a status attribute that the workflow uses to determine
        # which nodes are ready to be executed.
        data['status'] = 'new'

        self._backup_node(node_id, data)

        return self.graph.add_node(node_id, **data)

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

    def tasks(self, *args, **kwargs):
        """ Access tasks in the workflow. """
        return self.graph.nodes(*args, **kwargs)

    def get_task(self, task_id, fallback=None, data=False):
        """ Get a single task if it is present in the workflow.

        Similar to dict.get, this returns "fallback" if the task_id is not in
        the graph. """
        if task_id in self.graph:
            if data:
                return self.graph.nodes()[task_id]
            else:
                return task_id
        else:
            return fallback

    def send(self, task_id, return_code):
        if return_code is None:
            self.reset(task_id)
        elif return_code != 0:
            self.fail(task_id)
        else:
            self.complete(task_id)

    def complete(self, task_id, return_code=0):
        """ Complete a task.

        :param task_id: str task id
        :param return_code: int
        :return: None
        """
        log.critical('Task complete: {}'.format(task_id))

        self.update(
            task_id,
            return_code=return_code,
            datetime_end=str(datetime.now()),
            status='complete',
        )

    def fail(self, task_id, return_code=1):
        """ Fail a task. This also fails any tasks dependent on task_id.

        :param task_id: str task id
        :param return_code: int
        :return: None
        """
        log.critical('Task failed: {}'.format(task_id))

        self.update(
            task_id,
            return_code=return_code,
            datetime_end=str(datetime.now()),
            status='failed',
        )

        for d in self.graph.predecessors(task_id):
            self.fail(d, -42)

    @autosave
    def update(self, task_id, **kwargs):
        """ Change the status of a node. """
        self.last_update = utils.fingerprint()
        self.graph.nodes[task_id].update(**kwargs)

    @autosave
    def add_task(self, task_id, **kwargs):
        """ Add a task to the workflow. """
        self._add_node(task_id, kwargs)

    @autosave
    def add_dependency(self, task_id, before=None, after=None):
        """ Add a dependency to a task. """
        if task_id not in self.graph:
            raise ValueError('{} not in graph'.format(task_id))

        stack = []
        try:
            if before:
                if isinstance(before, str):
                    before = (before,)

                before_pats = [re.compile('^{}$'.format(b)) for b in before]

                for pat in before_pats:
                    matches = filter(pat.match, self.tasks())

                    if not matches:
                        raise ValueError('No matching tasks for: {}'.format(
                            pat.pattern))

                    for child_id in matches:
                        self._add_edge(from_node=child_id, to_node=task_id)
                        stack.append((child_id, task_id))

            if after:
                if isinstance(after, str):
                    after = (after,)

                after_pats = [re.compile('^{}$'.format(a)) for a in after]

                for pat in after_pats:
                    matches = filter(pat.match, self.tasks())

                    if not matches:
                        raise ValueError('No matching tasks for {}'.format(
                            pat.pattern))

                    for parent_id in matches:
                        self._add_edge(from_node=task_id, to_node=parent_id)
                        stack.append((task_id, parent_id))

        except Exception as e:
            for u, v in stack:
                log.debug('Rolling back: {} -> {}'.format(u, v))
                self.graph.remove_edge(u, v)
            raise e from None

    def task_ready(self, task_id):
        """ Returns True if the given "node_id" is ready for execution. """
        task_directives = self.graph.nodes()[task_id]

        if task_directives['status'] != 'new':
            return False

        for dependency in self.graph.successors(task_id):
            dependency_status = self.graph.nodes()[dependency]['status']
            if dependency_status != 'complete':
                return False
        else:
            log.critical('Task ready: {}'.format(task_id))
            return task_id, task_directives

    def status(self, task_id=None):
        """ Returns the status of a task. Returns status of all tasks if
        "task_id" is None. """
        if task_id is not None:
            return self.graph.nodes[task_id]['status']
        else:
            return [d['status'] for n, d in self.graph.nodes(data=True)]

    def compose(self, new_wf):
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
        existing_task_backup = self.graph.graph['_backup'].copy()
        new_task_backup = new_wf.graph.graph['_backup'].copy()
        existing_task_backup.update(new_task_backup)

        self.graph = nx.compose(new_wf.graph, self.graph)
        self.graph.graph['_backup'] = existing_task_backup

    def compose_all(self, *wfs):
        """ Compose this workflow with multiple other workflows.
        See Workflow.compose for more details. """
        wfs = [wf.graph for wf in wfs] + [self.graph]
        self.graph = nx.compose_all(wfs)

    def reset(self, task_id):
        """ Reset a task. This destroys any existing workflow records for the
        task. """
        log.critical('Task reset: {}'.format(task_id))
        self._restore_node(task_id)

    def reset_all(self):
        """ Reset all tasks. """
        for task_id, task_data in self.tasks(data=True):
            self.reset(task_id)

    def retry(self):
        """ Resets all pending and failed tasks. """
        for task_id, task_data in self.tasks(data=True):
            if task_data['status'] in ('pending', 'failed'):
                self.reset(task_id)

    def resume(self):
        """ Resets all pending tasks. """
        for task_id, task_data in self.tasks(data=True):
            if task_data['status'] == 'pending':
                self.reset(task_id)

    def save(self, path=None):
        if path is None:
            path = self.path or getattr(self.project, 'workflow_path', None)

            if path is None:
                raise ValueError('No path has been set for this workflow!')

        return save(self, path)


class WorkflowIterator(object):
    def __init__(self, workflow):
        self.workflow = workflow
        self._new_tasks = list()
        self._pending_tasks = list()
        self._complete_tasks = list()
        self._setup_lists()

    def __repr__(self):
        return '<WorkflowIterator {}>'.format(self.status())

    def _setup_lists(self):
        self._new_tasks = list()
        self._pending_tasks = list()
        self._complete_tasks = list()

        for task_id, task_data in self.workflow.tasks(True):
            status = task_data['status']

            if status == 'new':
                self._new_tasks.append(task_id)
            elif status == 'pending':
                self._pending_tasks.append(task_id)
            else:
                self._complete_tasks.append(task_id)

    def __next__(self):
        try:
            task_id = self._new_tasks.pop(0)
            task = self.workflow.task_ready(task_id)

            if task:
                self.workflow.update(task_id, status='pending')
                self._pending_tasks.append(task_id)
                return task
            else:
                self._new_tasks.append(task_id)

        except IndexError:
            if self._pending_tasks:
                return None
            else:
                raise StopIteration from None

    def __send__(self, task_id, returncode):
        self._pending_tasks.remove(task_id)
        res = self.workflow.send(task_id, returncode)

        if returncode != 0:
            # When a task fails, the workflow class is smart enough to fail
            # all dependencies of the task. But this iterator is not. So,
            # we need to request the new task status from the workflow after
            # failing a task.
            self._setup_lists()
        else:
            self._complete_tasks.append(task_id)

        return res

    def send(self, *args, **kwargs):
        return self.__send__(*args, **kwargs)

    def status(self):
        return {
            'new': len(self._new_tasks),
            'pending': len(self._pending_tasks),
            'complete': len(self._complete_tasks)
        }


def from_node_link_data(data):
    graph = json_graph.node_link_graph(data)
    return Workflow(graph=graph)


def to_node_link_data(wf):
    return json_graph.node_link_data(wf.graph)


def to_cytoscape_json_data(wf):
    """ Export a workflow as a cytoscape JSON file

    Cytoscape is good for vizualizing network graphs. It complains about
    node data that are not strings, so all node data are converted to string
    on export. This causes the cytoscape json files to be a one-way
    export, they cannot be loaded back into a workflow. """
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
    built. """
    log.critical('Building workflow...')

    if isinstance(tasks, str):
        # If tasks are not loaded yaml yet, do it automatically
        tasks = utils.yaml_loads(tasks)

    if not tasks:
        raise ValueError('No tasks were ')

    wf = Workflow()

    for task in tasks:
        wf.add_task(task_id=task['id'], **task)

    for task in tasks:
        # Before/After linking
        if 'before' in task:
            wf.add_dependency(task['id'], before=task['before'])

        if 'after' in task:
            wf.add_dependency(task['id'], after=task['after'])

        # In/Out linking
        # if 'in' in node:
        #     for dep in node['in']:
        #         found_match = False
        #         for n in wf.nodes(data=True):
        #             if n['out'] == dep:
        #                 wf.add_dependency(node['id'], after=(n,))
        #                 found_match = True
        #         if found_match is False:
        #             msg = 'Unmatched "in": {}'.format(node)
        #             raise WorkflowBuilderError(msg)

        # TODO This type of linking is much more complicated
        # need to revisit and think about how to address edge
        # cases: out directive with no ins, number of matches
        # allowed per in/out etc.
    
    log.critical('Workflow ready: {}'.format(wf))
    return wf
