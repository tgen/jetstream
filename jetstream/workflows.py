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
import logging
import pickle
import random
import re
import shutil
import textwrap
import traceback
from collections import Counter, deque
from datetime import datetime
from threading import Lock

import networkx as nx
from networkx.readwrite import json_graph

import jetstream
from jetstream import utils
from jetstream.tasks import Task


log = logging.getLogger(__name__)
workflow_extensions = {
    '': 'pickle',
    '.pickle': 'pickle',
    '.yaml': 'yaml',
    '.yml': 'yaml',
    '.json': 'json',
}


class NotDagError(ValueError):
    """Raised when edges are added that would result in a graph that is not
    a directed-acyclic graph"""


class Workflow(object):
    def __init__(self, tasks=None, **kwargs):
        built_w_version = kwargs.pop('jetstream_version', jetstream.__version__)
        current_version = jetstream.__version__

        if built_w_version != current_version:
            # TODO Only warn if built_w_version is higher than current?
            msg = f'This workflow was built with a different jetstream ' \
                  f'version:\n' \
                  f'current: {current_version}\n' \
                  f'workflow: {built_w_version}'
            log.warning(msg)

        self.graph = nx.DiGraph(jetstream_version=current_version, **kwargs)
        self.save_path = None
        self._lock = Lock()
        self._cm_stack = list()
        self._iter_tasks = list()
        self._iter_pending = list()

        if tasks:
            self.add_tasks(tasks)

    def __contains__(self, item):
        if isinstance(item, Task):
            return self.graph.__contains__(item.name)
        return self.graph.__contains__(item)

    def __enter__(self):
        """Workflows can be edited in a transaction using the context manager
        statement "with". This allows multiple task additions to take place
        with only a single update to the workflow edges. """
        self.lock()
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
        self.unlock()

    def __iter__(self):
        """In order to reduce the search time for the next available task,
        they are stored in separate lists, and then removed as they are
        completed. When a change to the graph occurs during iteration, these
        lists should be recalculated. """
        log.debug('Building workflow iterator...')
        self._iter_tasks = list()
        self._iter_pending = list()
        self._iter_done = list()

        for name in nx.topological_sort(self.graph):
            if self.get_task(name).is_new():
                self._iter_tasks.append(name)
            elif self.get_task(name).is_pending():
                self._iter_pending.append(name)
            else:
                self._iter_done.append(name)

        return self

    def __len__(self):
        return len(self.graph)

    def __next__(self):
        """Select the next available task for execution. If no task is ready,
        this will return None."""
        log.debug('Request for next task!')
        log.debug(f'{len(self._iter_tasks)} tasks remaining')
        log.debug(f'{len(self._iter_pending)} tasks pending')
        log.debug(f'{len(self._iter_done)} tasks done')

        if self.is_locked():
            raise RuntimeError('Workflow.__next__() called while locked!')

        # Drop all pending tasks that have completed since the last call
        _temp = list()
        for name in self._iter_pending:
            t = self.get_task(name)
            if t.is_done():
                self._iter_done.append(name)
            elif t.is_new():
                self._iter_tasks.append(name)
            else:
                _temp.append(name)
        self._iter_pending = _temp
        
        # Start search for next task
        for i in reversed(range(len(self._iter_tasks))):
            name = self._iter_tasks[i]
            task = self.get_task(name)

            if task.is_done():
                self._iter_tasks.pop(i)
                self._iter_done.append(name)
            elif task.is_pending():
                self._iter_tasks.pop(i)
                self._iter_pending.append(name)
            elif task.is_ready():
                self._iter_tasks.pop(i)
                self._iter_pending.append(name)
                task.pending()
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
        log.debug('Adding edge: {} -> {}'.format(from_node, to_node))

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

        if after:
            log.debug('Adding "after" edges for: {}'.format(task))
            matches = set()

            if isinstance(after, str):
                matches.add(after)
            elif isinstance(after, (list, tuple)):
                for target in after:
                    matches.add(target)
            elif isinstance(after, dict) and 're' in after:
                for target in self.find(after['re'], fallback=set()):
                    matches.add(target)
            else:
                raise ValueError(f'Unsupported "after" type in {task}')

            if task.name in matches:
                matches.remove(task.name)

            for match_name in matches:
                self._add_edge(from_node=match_name, to_node=task.name)

    def _make_edges_before(self, task):
        """Generate edges for "before" directives of a task
        "before" specifies edges that should run:

            task -------> tasks with name matching "before" pattern, ...

        """
        before = task.directives().get('before')

        if before:
            log.debug('Adding "before" edges for: {}'.format(task))
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

            if task.name in matches:
                matches.remove(task.name)

            for match_name in matches:
                self._add_edge(from_node=task.name, to_node=match_name)

    def _make_edges_input(self, task):
        """Generate edges for "input" directives of a task
        "input" specifies edges that should run:

            tasks with output matching "input" pattern, ... -------> task

        Where target includes an "output" value matching the "input" value."""
        input = task.directives().get('input')

        if input:
            log.debug('Adding "input" edges for: {}'.format(task))
            if isinstance(input, str):
                matches = self.find_by_output(input)
            elif isinstance(input, (list, tuple)):
                matches = set()
                for target in input:
                    for name in self.find_by_output(target):
                        matches.add(name)
            elif isinstance(input, dict) and 're' in input:
                # Inputs must be searched for in all nodes, but this allows
                # the syntax to match before/after directives.
                matches = self.find_by_output(input['re'])
            else:
                raise ValueError(f'Unsupported "input" type in {task}')

            if task.name in matches:
                matches.remove(task.name)

            for match_name in matches:
                self._add_edge(from_node=match_name, to_node=task.name)

    def _format_pattern(self, pat):
        """Pads a regex pattern so that it only matches exact strings"""
        return re.compile('^{}$'.format(pat))

    def _remove_node(self, name):
        """Directly remove nodes from the workflow graph. This can easily
        result in a corrupted workflow by leaving orphaned dependencies behind.
        Use Workflow.remove_tasks() instead."""
        log.debug(f'Removing node: {name}')
        self.graph.remove_node(name)

    def ancestors(self, task):
        """Returns a generator that yields all of the ancestors of a given task
        object or id. See also Workflow.dependencies() """
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.name

        return (self.get_task(name) for name in nx.ancestors(self.graph, task_id))

    def add_task(self, task=None, **kwargs):
        """Add a node to the graph and calculate any dependencies.

        Nodes are expected to be an instance of Jetstream.Task. If the workflow
        is not locked (via "with" statement) this will trigger Workflow.update.
        """
        if task is None:
            task = Task(**kwargs)

        if task.name in self.graph:
            raise ValueError('Duplicate task ID: {}'.format(task.name))

        log.debug('Adding task: {}'.format(task))

        if task.workflow:
            task = task.copy()

        task.workflow = self
        self.graph.add_node(task.name, obj=task)

        if self.is_locked():
            self._cm_stack.append(task.name)
        else:
            try:
                self.update()
            except Exception as e:
                self.graph.remove_node(task.name)
                raise e

        return task

    def add_tasks(self, tasks):
        with self:
            for task in tasks:
                try:
                    if isinstance(task, Task):
                        self.add_task(task)
                    else:
                        self.add_task(**task)
                except Exception:
                    tb = traceback.format_exc()
                    txt = textwrap.shorten(str(task), 42)
                    msg = f'Error with loading task: {txt}\n\n' \
                          f'{textwrap.indent(tb, "  ")}\n'
                    raise ValueError(msg) from None

    def complete_tasks(self, pattern, *args, **kwargs):
        log.debug(f'Complete tasks: {pattern}')

        for task in self.find(pattern, objs=True):
            task.complete(*args, **kwargs)

    def dependencies(self, task):
        """Returns a generator that yields all only the direct dependencies of
        a given task object or id. See also Workflow.ancestors() """
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.name

        return (self.get_task(name) for name in self.graph.predecessors(task_id))

    def dependents(self, task):
        """Returns a generator that yields only the tasks that directly
        depend upon of a given task object or id. See also
        Workflow.decendants() """
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.name

        return (self.get_task(name) for name in self.graph.successors(task_id))

    def descendants(self, task):
        """Returns a generator that yields all the descendants of a given
        task object or id. See also Workflow.dependents() """
        if isinstance(task, str):
            task_id = task
        else:
            task_id = task.name

        return (self.get_task(name) for name in nx.descendants(self.graph, task_id))

    def draw(self, *args, **kwargs):
        """Attempt to draw the network graph for this workflow. See
        jetstream.workflows.draw_workflow() for more info."""
        return draw_workflow(self, *args, **kwargs)

    def find(self, pattern, fallback=utils.sentinel, objs=False):
        """Find tasks by matching the pattern with the task ID. If no matches
        are found, and fallback is not set, a ValueError will be raised. If
        fallback is set, it will be returned when no matches are found."""
        log.debug('Find: {}'.format(pattern))

        pat = self._format_pattern(pattern)
        matches = set([name for name in self.graph.nodes() if pat.match(name)])

        if matches:
            if objs:
                return set([self.get_task(t) for t in matches])
            else:
                return matches
        elif fallback is utils.sentinel:
            if pattern == '.*':
                return set()
            raise ValueError('No task names match value: {}'.format(pattern))
        else:
            return fallback

    def find_by_output(self, pattern, fallback=utils.sentinel, objs=False):
        """Find tasks where the pattern matches the task output directives. If
        no matches are found, and fallback is not set, a ValueError will be
        raised. If fallback is set, it will be returned when no matches are
        found."""
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
            if objs:
                return set([self.get_task(t) for t in matches])
            else:
                return matches
        elif fallback is utils.sentinel:
            raise ValueError('No task outputs match pattern: {}'.format(pattern))
        else:
            return fallback

    def get_task(self, task_id):
        """Return the task object for a given task_id"""
        return self.graph.nodes[task_id]['obj']

    def lock(self):
        self._lock.acquire()

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
        """Returns all task objects as a list"""
        return list(self.tasks(objs=True))

    def mash(self, workflow):
        """Mash this workflow with another."""
        return mash(self, workflow)

    def new_task(self,  **kwargs):
        """Shortcut to create a new Task object and add to this workflow."""
        task = Task(**kwargs)
        return self.add_task(task)

    def remove_tasks(self, pattern, *, descendants=False, force=False):
        """Remove task(s) from the workflow.
        This will find tasks by name and remove each match. """
        log.debug(f'Remove tasks: {pattern}')

        for task in self.find(pattern, objs=True):
            try:
                has_deps = next(task.dependents())
            except StopIteration:
                has_deps = False

            if has_deps:
                if descendants:
                    for d in task.descendants():
                        self._remove_node(d.name)
                    self._remove_node(task.name)
                elif force:
                    self._remove_node(task.name)
                else:
                    err = f'{task} has decendants that would be orphaned if ' \
                          'it were removed. Use "descendants=True" to remove' \
                          'this task and any descendants.'
                    raise ValueError(err)
            else:
                self._remove_node(task.name)

    def reset_tasks(self, pattern, *args, **kwargs):
        log.debug(f'Reset tasks: {pattern}')

        for task in self.find(pattern, objs=True):
            task.reset(*args, **kwargs)

    def fail_tasks(self, pattern, *args, **kwargs):
        log.debug(f'Fail tasks: {pattern}')

        for task in self.find(pattern, objs=True):
            task.reset(*args, **kwargs)

    def reset(self, method):
        """Resets state for tasks in this workflow
        Resetting tasks allows workflows that were stopped prior to completion
        to be rerun. Or tasks that failed due to external state can be retried:

        all - Resets state for all tasks
        resume - Resets state for all "pending" tasks
        retry - Resets state for all "pending" and "failed" tasks
        """
        if method == 'retry':
            self.retry()
        elif method == 'resume':
            self.resume()
        elif method == 'all':
            self.reset_all()
        else:
            err = f'Unrecognized workflow reset method: {method} Choose one ' \
                  f'of all, retry, resume'
            raise ValueError(err)

    def reset_all(self):
        """Resets state for all tasks"""
        log.critical('Reset: Resetting state for all tasks...')
        for task in self.tasks(objs=True):
            task.reset(descendants=False)

    def resume(self):
        """Resets state for all "pending" tasks """
        log.info('Resume: Resetting state for all pending tasks...')
        for task in self.tasks(objs=True):
            if task.status == 'pending':
                task.reset(descendants=False)

    def retry(self):
        """Resets state for all "pending" and "failed" tasks """
        log.info('Retry: Resetting state for all pending and failed tasks...')
        for task in self.tasks(objs=True):
            if task.status in ('pending', 'failed'):
                task.reset(descendants=False)

    def to_dict(self):
        """Convert the workflow to a node-link formatted object that can
        be easily dumped to JSON/YAML """
        return to_dict(self)

    @staticmethod
    def from_dict(data):
        """Given node-link dictionary, generate a workflow object"""
        return from_dict(data)

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
        If objs is False, only the names will be returned."""
        if objs:
            return (t['obj'] for i, t in self.graph.nodes(data=True))
        else:
            return self.graph.nodes()

    def to_json(self, indent=None, sort_keys=True):
        """Returns JSON string representation of this workflow"""
        log.debug('Dumping to json...')
        return utils.json_dumps(
            self.to_dict(), indent=indent, sort_keys=sort_keys
        )

    def to_yaml(self):
        """Returns YAML string representation of this workflow"""
        s = self.to_dict()
        log.debug('Dumping to yaml...')
        return utils.yaml_dumps(s)

    def update(self):
        """Recalculate the edges for this workflow"""
        log.debug('Adding edges...')
        for task in self.tasks(objs=True):
            self._make_edges_after(task)
            self._make_edges_before(task)
            self._make_edges_input(task)

    def unlock(self):
        self._lock.release()


def to_dict(workflow):
    log.debug('Converting to node link data...')
    data = json_graph.node_link_data(workflow.graph)

    log.debug('Serializing nodes...')
    for node in data['nodes']:
        node['obj'] = node['obj'].to_dict()
        node['obj'].pop('name')

    return data


def from_dict(data):
    wf = Workflow(**data['graph'])

    with wf:
        for node in data['nodes']:
            task_data = node['obj']
            task_data['name'] = node['id']
            t = jetstream.tasks.from_dict(task_data)
            wf.add_task(t)

    return wf


def draw_workflow(wf, *, figsize=(12, 12), cm=None, filename=None,
                  interactive=False, **kwargs):
    """This requires matplotlib setup and graphviz installed. It is not
    very simple to get these dependencies setup, so they're loaded as needed
    and not required for package install."""
    try:
        import matplotlib
        if not interactive and os.environ.get('DISPLAY', '') == '':
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

    # labels = {id: node['obj'].label for id, node in wf.graph.nodes(data=True)}
    colors = [cm[node['obj'].status] for id, node in
              wf.graph.nodes(data=True)]

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


# TODO Move these to the config file
def get_workflow_loaders():
    return {
        'pickle': load_workflow_pickle,
        'yaml': load_workflow_yaml,
        'json': load_workflow_json
    }


def get_workflow_savers():
    return {
        'pickle': save_workflow_pickle,
        'yaml': save_workflow_yaml,
        'json': save_workflow_json
    }


def load_workflow(path, format=None):
    """Load a workflow from a file.

    This helper function will try to choose the correct file format based
    on the extension of the path, but defaults to pickle for unrecognized
    extensions. It also sets workflow.save_path to the path"""
    if format is None:
        ext = os.path.splitext(path)[1]
        format = workflow_extensions.get(ext, 'pickle')

    loader_fn = get_workflow_loaders()[format]
    log.debug(f'Loading workflow from: {path} with: {loader_fn}')

    wf = loader_fn(path)
    wf.save_path = os.path.abspath(path)
    return wf


def load_workflow_yaml(path):
    data = utils.load_yaml(path)
    return Workflow.from_dict(data)


def load_workflow_json(path):
    data = utils.load_json(path)
    return Workflow.from_dict(data)


def load_workflow_pickle(path):
    with open(path, 'rb') as fp:
        data = pickle.load(fp)
    return Workflow.from_dict(data)


def mash(G, H):
    """Mash together two Workflows

    Description::

        Example 1: Tasks already present remain completed
                   G         --->      H        =  New Workflow
          ---------------------------------------------------
                                   task1(new)       task1(complete)
            task1(complete)  --->      |        =     |
                                   task2(new)       task2(new)


        Example 2: Descendants of new or modified tasks are reset
                    G        --->      H       =  New Workflow
          ---------------------------------------------------
                                   task1(new)       task1(new)
            task2(complete)  --->      |        =     |
                                   task2(new)       task2(new)

    :param G: Workflow to be mashed, workflow properties will be kept from G
    :param H: Another Workflow
    :return: Workflow
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
                existing_task = wf.get_task(task.name)
            except KeyError:
                t = wf.add_task(task)
                new.add(t)
                continue

            # For tasks in H that are also in G replace if they've been modified
            if task.identity != existing_task.identity:
                wf.remove_tasks(existing_task.name, force=True)
                t = wf.add_task(task)
                modified.add(t)

    log.debug('Identifying tasks that need to be reset...')
    new_u_modified = new.union(modified)
    to_reset = set()
    temp_graph = wf.graph.copy()

    for t in new_u_modified:
        if t.name not in temp_graph:
            continue

        for d in nx.descendants(temp_graph, t.name):
            to_reset.add(d)
            temp_graph.remove_node(d)

    for t in to_reset:
        wf.get_task(t).reset(descendants=False)

    log.info(
        'Mash report:\n'
        f'New tasks: {len(new)}\n'
        f'Modified tasks: {len(modified)}\n'
        f'Reset due to modified ancestor: {len(to_reset)}'
    )

    return wf


def random_workflow(n=50, timeout=None, connectedness=3, trail=10):
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
    added = 0
    start = datetime.now()
    task = jetstream.tasks.random_task()
    queue = deque(maxlen=trail)
    wf.add_task(task)
    queue.append(task)

    while 1:
        if n and added >= n:
            log.critical('Task limit reached!')
            break

        if timeout and (datetime.now() - start).seconds > timeout:
            log.critical('Timeout reached!')
            break

        conns = random.randint(0, connectedness)
        inputs = []
        for i in range(conns):
            task = random.choice(queue)
            output = task.directives().get('output')
            inputs.append(output)

        try:
            task = jetstream.tasks.random_task(input=inputs)
            wf.add_task(task)
            queue.append(task)
            log.info('{} tasks added!'.format(added))
            added += 1
        except Exception as e:
            log.exception(e)

    return wf


def save_workflow(workflow, path, format=None):
    """Save a workflow to the path

    This helper function will try to choose the correct file format based
    on the extension of the path, but defaults to pickle for unrecognized
    extensions.

    :param workflow: Workflow instance
    :param path: where to save
    :param format: format to save
    :return: None
    """
    if format is None:
        ext = os.path.splitext(path)[1]
        format = workflow_extensions.get(ext, 'pickle')

    start = datetime.now()
    get_workflow_savers()[format](workflow, path)
    elapsed = datetime.now() - start

    log.debug('Workflow saved (after {}): {}'.format(elapsed, path))


def save_workflow_yaml(workflow, path):
    log.debug('Saving workflow (yaml): {}'.format(path))
    lock_path = path + '.lock'

    with open(lock_path, 'w') as fp:
        utils.yaml_dump(workflow.to_dict(), fp)

    shutil.move(lock_path, path)


def save_workflow_json(workflow, path):
    log.debug('Saving workflow (json): {}'.format(path))
    lock_path = path + '.lock'

    with open(lock_path, 'w') as fp:
        utils.json_dump(workflow.to_dict(), fp)

    shutil.move(lock_path, path)


def save_workflow_pickle(workflow, path):
    log.debug('Saving workflow (pickle): {}'.format(path))
    lock_path = path + '.lock'

    with open(lock_path, 'wb') as fp:
        pickle.dump(workflow.to_dict(), fp)

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

