"""Directed-acyclic graph models of computational workflows

The `Workflow` class models computational workflows as a directed-acyclic graph,
where nodes are tasks to complete, and edges represent dependencies between
those tasks. It includes methods for building workflows (new_task, add_task,
remove_task) in addition to the methods required for executing a workflow
(__next__, __iter__, dependencies, etc.).

A unique feature of Jetstream workflows is that task dependencies are
recalculated after every change to the workflow. Task dependencies (edges) are
not set up individually by the user, they're calculated from flow directives
that can match one or many tasks in the workflow.

Workflows can also be built by rendering a template. Templates are text
documents that describe a set of tasks to complete. They can include complex
dynamic elements via the templating syntax. Data used to render the templates
can be saved in files located in the config directory of project, or given as
arguments to template.

"""
import os
import re
import shutil
import pickle
import random
from datetime import datetime
import networkx as nx
from networkx.readwrite import json_graph
from threading import Lock
from collections import Counter
from pkg_resources import get_distribution
import jetstream
from jetstream import utils, log
from jetstream.tasks import Task

__version__ = get_distribution('jetstream').version


class NotDagError(ValueError):
    """Raised when edges are added that would result in a graph that is not
    a directed-acyclic graph"""

class Workflow(object):
    def __init__(self, **kwargs):
        built_w_version = kwargs.pop('jetstream_version', __version__)

        if built_w_version != __version__:
            # TODO Only warn if built_w_version is higher than current
            log.warning('This workflow was built with a different version')

        self.graph = nx.DiGraph(jetstream_version=__version__, **kwargs)
        self.save_path = None
        self._lock = Lock()
        self._cm_stack = list()
        self._iter_tasks = list()
        self._iter_pending = list()

    def __contains__(self, item):
        if isinstance(item, Task):
            return self.graph.__contains__(item.tid)
        return self.graph.__contains__(item)

    def __enter__(self):
        """Workflows can be edited in a transaction using the context manager
        statement "with". This allows multiple task additions to take place
        with only a single update to the workflow edges. """
        self._lock.acquire()
        self._cm_stack = list()

    def __exit__(self, exc_type, exc_value, traceback):
        """If there is an error during the transaction, all nodes will
        be rolled back on exit. """
        if exc_value is not None:
            for task_id in self._cm_stack:
                self.graph.remove_node(task_id)

        self.update()
        if not nx.is_directed_acyclic_graph(self.graph):
            for task_id in self._cm_stack:
                self.graph.remove_node(task_id)
            raise NotDagError

        self._cm_stack = list()
        self._lock.release()

    def __iter__(self):
        """In order to reduce the search time for the next available task,
        they are stored in separate lists, and then removed as they are
        completed. When a change to the graph occurs during iteration, these
        lists should be recalculated. """
        log.verbose('Building workflow iterator...')

        self._iter_tasks = list()
        self._iter_pending = list()
        self._iter_done = list()

        for task in nx.topological_sort(self.graph):
            if self.get_task(task).is_new():
                self._iter_tasks.append(task)
            elif self.get_task(task).is_pending():
                self._iter_pending.append(task)
            else:
                self._iter_done.append(task)

        return self

    def __len__(self):
        return len(self.graph)

    def __next__(self):
        """Select the next available task for execution. If no task is ready,
        this will return None."""
        log.verbose('Request for next task!')
        log.verbose(f'{len(self._iter_tasks)} tasks remaining')
        log.verbose(f'{len(self._iter_pending)} tasks pending')
        log.verbose(f'{len(self._iter_done)} tasks done')

        if self.is_locked():
            raise RuntimeError('Workflow.__next__() called while locked!')

        # Drop all pending tasks that have completed since the last call
        self._iter_pending = [t for t in self._iter_pending if not t.is_done()]
        
        # Search for the next available task 
        for i in reversed(range(len(self._iter_tasks))):
            tid = self._iter_tasks[i]
            task = self.get_task(tid)

            if task.is_done():
                self._iter_tasks.pop(i)
                self._iter_done.append(task)
            elif task.is_pending():
                self._iter_tasks.pop(i)
                self._iter_pending.append(task)
            elif task.is_ready():
                self._iter_tasks.pop(i)
                self._iter_pending.append(task)
                task.pending(quiet=True)
                return task
        
        # If there are any remaining or pending, return None until one is ready
        if self._iter_tasks or self._iter_pending:
            return None
        else:
            raise StopIteration

    def __repr__(self):
        stats = Counter([t.status for t in self.tasks(objs=True)])
        return '<jetstream.Workflow {}>'.format(stats)

    def _add_edge(self, from_node, to_node):
        """Edges represent dependencies between tasks. Edges run FROM one node
        TO another dependent node. Nodes can have multiple edges, but not
        multiple instances of the same edge (multigraph).

            Parent ---- Is a dependency of --->  Child
         (from_node)                           (to_node)

        This means that the in-degree of a node represents the number of
        dependencies it has. A node with zero in-edges is a "root" node, or a
        task with no dependencies. """
        log.verbose('Adding edge: {} -> {}'.format(from_node, to_node))

        if from_node not in self.graph:
            raise NotDagError(f'{from_node} is not in the workflow!')

        if to_node not in self.graph:
            raise NotDagError(f'{to_node} is not in the workflow!')

        self.graph.add_edge(from_node, to_node)

        if not self.is_locked():
            if not nx.is_directed_acyclic_graph(self.graph):
                self.graph.remove_edge(from_node, to_node)
                raise NotDagError('{} -> {}'.format(from_node, to_node))

    def _make_edges_after(self, task):
        """Generate edges for "after" directives of a task.
        "after" directives create edges that run:

            tasks with name matching "after" pattern, ...  ------>  task

        """
        after = task.directives().get('after')
        log.verbose('"after" directive: {}'.format(after))

        if after:
            matches = set()

            if isinstance(after, str):
                matches.add(after)
            elif isinstance(after, (list, tuple)):
                for target in after:
                    matches.add(target)
            elif isinstance(after, dict) and 're' in after:
                for target in self.find(after['re']):
                    matches.add(target)
            else:
                raise ValueError(f'Unsupported "after" type in {task}')

            if task.tid in matches:
                matches.remove(task.tid)

            for match_tid in matches:
                self._add_edge(from_node=match_tid, to_node=task.tid)

    def _make_edges_before(self, task):
        """Generate edges for "before" directives of a task
        "before" specifies edges that should run:

            task -------> tasks with name matching "before" pattern, ...

        """
        before = task.directives().get('before')
        log.verbose('"before" directive: {}'.format(before))

        if before:
            matches = set()

            if isinstance(before, str):
                matches.add(before)
            elif isinstance(before, (list, tuple)):
                for target in before:
                    matches.add(target)
            elif isinstance(before, dict) and 're' in before:
                for target in self.find(before['re']):
                    matches.add(target)
            else:
                raise ValueError(f'Unsupported "before" type in {task}')

            if task.tid in matches:
                matches.remove(task.tid)

            for match_tid in matches:
                self._add_edge(from_node=task.tid, to_node=match_tid)

    def _make_edges_input(self, task):
        """Generate edges for "input" directives of a task
        "input" specifies edges that should run:

            tasks with output matching "input" pattern, ... -------> task

        Where target includes an "output" value matching the "input" value."""
        input = task.directives().get('input')
        log.verbose('"input" directive: {}'.format(input))

        if input:
            if isinstance(input, str):
                matches = self.find_by_output(input)
            elif isinstance(input, (list, tuple)):
                matches = set()
                for target in input:
                    for tid in self.find_by_output(target):
                        matches.add(tid)
            elif isinstance(input, dict) and 're' in input:
                # Inputs must be searched for in all nodes, but this allows
                # the syntax to match before/after directives.
                matches = self.find_by_output(input['re'])
            else:
                raise ValueError(f'Unsupported "input" type in {task}')

            if task.tid in matches:
                matches.remove(task.tid)

            for match_tid in matches:
                self._add_edge(from_node=match_tid, to_node=task.tid)

    def _format_pattern(self, pat):
        return re.compile('^{}$'.format(pat))

    def ancestors(self, task):
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.tid

        return (self.get_task(tid) for tid in nx.ancestors(self.graph, task_id))

    def add_task(self, task):
        """Add a node to the graph and calculate any dependencies.

        Nodes are expected to be an instance of Jetstream.Task. If the workflow
        is not locked (via "with" statement) this will trigger Workflow.update.
        """
        if not isinstance(task, Task):
            raise ValueError('task must be instance of {}'.format(Task))

        if task.tid in self.graph:
            raise ValueError('Duplicate task ID: {}'.format(task.tid))

        log.verbose('Adding task: {}'.format(task))

        if task.workflow:
            task = task.copy()

        task.workflow = self
        self.graph.add_node(task.tid, obj=task)

        if self.is_locked():
            self._cm_stack.append(task.tid)
        else:
            try:
                self.update()
            except Exception as e:
                self.graph.remove_node(task.tid)
                raise e

        return task

    def dependencies(self, task):
        """Returns a generator that yields the dependencies of a given task"""
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.tid

        return (self.get_task(tid) for tid in self.graph.predecessors(task_id))

    def dependents(self, task):
        """Returns a generator that yields the dependents of a given task"""
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.tid

        return (self.get_task(tid) for tid in self.graph.successors(task_id))

    def descendants(self, task):
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.tid

        return (self.get_task(tid) for tid in nx.descendants(self.graph, task_id))

    def draw(self, *args, **kwargs):
        return draw_workflow(self, *args, **kwargs)

    def find(self, pattern, fallback=utils.sentinel):
        log.debug('Find: {}'.format(pattern))

        pat = self._format_pattern(pattern)
        matches = set([tid for tid in self.graph.nodes() if pat.match(tid)])

        if matches:
            return matches

        if fallback is utils.sentinel:
            if pattern == '.*':
                return set()
            raise ValueError('No task names match value: {}'.format(pattern))
        else:
            return fallback

    def find_by_output(self, pattern, fallback=utils.sentinel):
        pat = self._format_pattern(pattern)
        matches = set()

        for task_id, data in self.graph.nodes(True):
            task = data['obj']

            try:
                output = task.directives()['output']
            except KeyError:
                continue

            output = utils.coerce_sequence(output)

            for value in output:
                if pat.match(value):
                    matches.add(task_id)
                    break

        if matches:
            return matches
        elif fallback is utils.sentinel:
            raise ValueError('No task outputs match pattern: {}'.format(pattern))
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

    def mash(self, workflow):
        """Mash this workflow with another."""
        return mash(self, workflow)

    def new_task(self, *args, **kwargs):
        """Shortcut to create a new Task object and add to this workflow."""
        task = Task(*args, **kwargs)
        return self.add_task(task)

    def remove_task(self, pattern, force=False, descendants=False):
        """Remove task(s) from the workflow.
        This will find tasks by name and call remove_task_id for each match. """
        log.info('Remove task: {}'.format(pattern))
        matches = self.find(pattern)

        for task_id in matches:
            try:
                deps = next(self.dependents(task_id))
            except StopIteration:
                deps = None

            if deps is None or force:
                self.graph.remove_node(task_id)
            else:
                raise ValueError('Task has dependents!')

    def reset(self):
        """Resets all tasks state."""
        log.critical('Resetting state for all tasks...')
        for task in self.tasks(objs=True):
            task.reset()

    def resume(self):
        """Resets all "pending" tasks state."""
        log.info('Resetting state for all pending tasks...')
        for task in self.tasks(objs=True):
            if task.status == 'pending':
                task.reset()

    def retry(self):
        """Resets all "pending" and "failed" tasks state."""
        log.info('Resetting state for all pending and failed tasks...')
        for task in self.tasks(objs=True):
            if task.status in ('pending', 'failed'):
                task.reset()

    def serialize(self):
        """Convert the workflow to a node-link formatted object that can
        be easily dumped to JSON/YAML """
        log.debug('Converting to node link data...')
        data = json_graph.node_link_data(self.graph)

        log.debug('Serializing nodes...')
        for node in data['nodes']:
            node['obj'] = node['obj'].serialize()
            node['obj'].pop('tid')

        return data

    @staticmethod
    def deserialize(data):
        """Given node-link data, generate a workflow object"""
        wf = Workflow(**data['graph'])

        with wf:
            for node in data['nodes']:
                task_data = node['obj']
                task_data['tid'] = node['id']
                wf.new_task(from_data=task_data)

        return wf

    def save(self, path=None, *args, **kwargs):
        """Shortcut to :func:`jetstream.save_workflow`"""
        path = path or self.save_path

        if not path:
            raise ValueError('No save path has been set')

        return save_workflow(self, path, *args, **kwargs)

    @staticmethod
    def load(*args, **kwargs):
        """Shortcut to :func:`jetstream.load_workflow`"""
        return load_workflow(*args, **kwargs)

    def tasks(self, objs=True):
        """Access the tasks in this workflow.
        If objs is False, only the tids will be returned."""
        if objs:
            return (t['obj'] for i, t in self.graph.nodes(data=True))
        else:
            return self.graph.nodes()

    def to_json(self, indent=None, sort_keys=True):
        """Returns JSON string representation of this workflow"""
        log.debug('Dumping to json...')
        return utils.json_dumps(
            self.serialize(), indent=indent, sort_keys=sort_keys
        )

    def to_yaml(self):
        """Returns YAML string representation of this workflow"""
        s = self.serialize()
        log.debug('Dumping to yaml...')
        return utils.yaml_dumps(s)

    def update(self):
        """Recalculate the edges for this workflow"""
        log.info('Updating workflow DAG...')
        for task in self.tasks(objs=True):
            self._make_edges_after(task)
            self._make_edges_before(task)
            self._make_edges_input(task)


def mash(G, H):
    """Mash together two workflows
    ::

    Example 1: Tasks already present remain completed
               G         --->      H        =     G
      ---------------------------------------------------
                               task1(new)       task1(complete)
        task1(complete)  --->      |        =     |
                               task2(new)       task2(new)


    Example 2: Descendants of new or modified tasks are reset
                G         --->      H        =     G
      ---------------------------------------------------
                               task1(new)       task1(new)
        task2(complete)  --->      |        =     |
                               task2(new)       task2(new)

    :param G: Workflow
    :param H: Workflow
    :return: None
    """
    log.info(f'Mashing {G} with {H}')

    wf = jetstream.Workflow(**G.graph.graph)
    new = set()
    modified = set()

    with wf:
        log.debug('Adding all tasks from G...')
        for task in G.tasks(objs=True):
            wf.add_task(task)

        log.debug('Adding tasks from H...')
        for task in H.tasks(objs=True):
            # Identify tasks in H that are not in G yet
            try:
                existing_task = wf.get_task(task.tid)
            except KeyError:
                t = wf.add_task(task)
                new.add(t)
                continue

            # For tasks in H that are also in G replace if they've been modified
            if task.identity != existing_task.identity:
                wf.remove_task(existing_task.tid, force=True)
                t = wf.add_task(task)
                modified.add(t)

    log.debug('Identifying tasks that need to be reset...')
    new_u_modified = new.union(modified)
    to_reset = set()
    temp_graph = wf.graph.copy()

    for t in new_u_modified:
        if t.tid not in temp_graph:
            continue

        for d in nx.descendants(temp_graph, t.tid):
            to_reset.add(d)
            temp_graph.remove_node(d)

    for t in to_reset:
        wf.get_task(t).reset(quiet=True, descendants=False)

    log.info(
        'Mash report:\n'
        f'New tasks: {len(new)}\n'
        f'Modified tasks: {len(modified)}\n'
        f'Ancestors indirectly affected: {len(to_reset)}'
    )

    return wf


def draw_workflow(wf, *, figsize=(12,12), cm=None, filename=None, **kwargs):
    """This requires matplotlib setup and graphviz installed. It is not
    very simple to get these dependencies setup, so they're loaded as needed
    and not required for package install."""
    try:
        import matplotlib
        if os.environ.get('DISPLAY','') == '':
            log.critical('No display found. Using non-interactive Agg backend')
            matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from networkx.drawing.nx_agraph import graphviz_layout
    except ImportError:
        log.critical('This feature requires matplotlib and graphviz')
        raise
    
    f = plt.figure(figsize=figsize)

    cm = cm or {
        'new': 'white',
        'pending': 'yellow',
        'failed': 'red',
        'complete': 'green'
    }

    #labels = {id: node['obj'].label for id, node in wf.graph.nodes(data=True)}
    colors = [cm[node['obj'].status] for id, node in wf.graph.nodes(data=True)]

    nx.draw(
        wf.graph,
        pos=graphviz_layout(wf.graph, prog='dot'),
        ax=f.add_subplot(111),
       # labels=labels,
        node_color=colors,
        linewidths=1,
        edgecolors='black',
        with_labels=True,
        **kwargs
    )

    if filename is not None:
        f.savefig(filename)

    return plt


def load_workflow(path, format=None):
    """Load a workflow from a file.

    This helper function will try to choose the correct file format based
    on the extension of the path, but defaults to pickle for unrecognized
    extensions. It also sets workflow.save_path to the path"""
    if format is None:
        ext = os.path.splitext(path)[1]
        format = jetstream.workflow_extensions.get(ext, 'pickle')

    wf = jetstream.workflow_loaders[format](path)
    wf.save_path = os.path.abspath(path)
    return wf


def load_workflow_yaml(path):
    data = utils.load_yaml(path)
    return Workflow.deserialize(data)


def load_workflow_json(path):
    data = utils.load_json(path)
    return Workflow.deserialize(data)


def load_workflow_pickle(path):
    with open(path, 'rb') as fp:
        data = pickle.load(fp)
    return Workflow.deserialize(data)


def save_workflow(workflow, path, format=None):
    """Save a workflow to the path

    This helper function will try to choose the correct file format based
    on the extension of the path, but defaults to pickle for unrecognized
    extensions.

    :param workflow: Workflow instance
    :param path: where to save
    :return: None
    """
    if format is None:
        ext = os.path.splitext(path)[1]
        format = jetstream.workflow_extensions.get(ext, 'pickle')

    start = datetime.now()
    jetstream.workflow_savers[format](workflow, path)
    elapsed = datetime.now() - start

    log.debug('Workflow saved (after {}): {}'.format(elapsed, path))


def random_workflow(n=50, timeout=None, connectedness=3, start_numbering=0):
    """Random workflow generator. The time to generate a random task for a
     workflow scales exponentially, so this can take a very long time for
     large numbers of tasks.

     Workflow size can be controlled by n, timeout, or both. But, at least
     one must be set. Timeout is the number of seconds (approx.) before the
     workflow will be returned.

     Connectedness is the maximum number of connections to proc when
     generating each task.

     """
    if not (n or timeout):
        raise ValueError('Must set n or timeout')

    wf = Workflow()
    cmds = (None, 'echo', 'hostname', 'ls', 'date', 'who', 'sleep 1')
    start = datetime.now()
    added = start_numbering

    if n:
        n = start_numbering + n

    while 1:
        if n and added >= n:
            log.critical('Task limit reached!')
            break

        if timeout and (datetime.now() - start).seconds > timeout:
            log.critical('Timeout reached!')
            break

        tasks = wf.list_tasks()
        directives = {}
        directives['name'] = str(added)
        directives['cmd'] = random.choice(cmds)
        directives['output'] = random.choice([
            None,
            hex(random.getrandbits(32)) + '.txt'
        ])

        if tasks:
            for i in range(random.randint(0, connectedness)):
                task = random.choice(tasks)
                name = task.directives().get('name')
                output = task.directives().get('output')

                if random.random() > 0.5:
                    if output is not None:
                        directives['input'] = output
                else:
                    directives['after'] = name

        try:
            wf.new_task(**directives)
            log.info('{} tasks added!'.format(added))
            added += 1
        except Exception as e:
            log.exception(e)

    return wf


def save_workflow_yaml(workflow, path):
    log.info('Saving workflow (yaml): {}'.format(path))
    lock_path = path + '.lock'

    with open(lock_path, 'w') as fp:
        utils.yaml_dump(workflow.serialize(), fp)

    shutil.move(lock_path, path)


def save_workflow_json(workflow, path):
    log.info('Saving workflow (json): {}'.format(path))
    lock_path = path + '.lock'

    with open(lock_path, 'w') as fp:
        utils.json_dump(workflow.serialize(), fp)

    shutil.move(lock_path, path)


def save_workflow_pickle(workflow, path):
    log.info('Saving workflow (pickle): {}'.format(path))
    lock_path = path + '.lock'

    with open(lock_path, 'wb') as fp:
        pickle.dump(workflow.serialize(), fp)

    shutil.move(lock_path, path)


def to_cytoscape_json_data(wf):
    """Export a workflow as a cytoscape JSON file

    Cytoscape is good for vizualizing network graphs. It complains about node
    data that are not strings, so all node data are converted to strings on
    export. This causes Cytoscape json files to be a one-way export, they
    cannot be loaded back to workflow objects. """
    data = nx.cytoscape_data(wf.graph)

    for n in data['elements']['nodes']:
        for k, v in n['data'].items():
            n['data'][k] = str(v)

    return data

