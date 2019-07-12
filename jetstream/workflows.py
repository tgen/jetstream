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
import logging
import glob
import pickle
import random
import re
import shutil
from collections import Counter, deque
from datetime import datetime
from distutils.version import LooseVersion
import networkx as nx
import jetstream
from jetstream import utils
from jetstream.tasks import Task

log = logging.getLogger(__name__)


class Workflow:
    def __init__(self, tasks=None, props=None, path=None, version=None):
        if tasks:
            self.tasks = {task.name: task for task in tasks}
        else:
            self.tasks = {}
        self.props = props or {}
        self.path = path
        self.version = version or jetstream.__version__

    def __contains__(self, item):
        if isinstance(item, str):
            return item in self.tasks
        else:
            return item.name in self.tasks

    def __getitem__(self, item):
        if isinstance(item, str):
            return self.tasks.__getitem__(item)
        else:
            return self.tasks.__getitem__(item.name)

    def __getstate__(self):
        """Enables checking versions when workflows are loaded"""
        return self.__dict__

    def __iter__(self):
        return iter(self.tasks.values())

    def __len__(self):
        return len(self.tasks)

    def __setstate__(self, d):
        """Enables checking versions when workflows are loaded"""
        self.__dict__ = d
        self.check_versions()

    def add(self, task):
        if task.name in self.tasks:
            err = f'Duplicate task name added to workflow: {task.name}'
            raise ValueError(err)
        self.tasks[task.name] = task

    def check_versions(self):
        current_version = LooseVersion(jetstream.__version__)
        built_w_version = LooseVersion(self.version)

        if built_w_version > current_version:
            msg = f'This workflow was built with a newer version of ' \
                  f'Jetstream: {built_w_version}'
            log.warning(msg)
            return False
        else:
            return True

    def find(self, pattern, style='regex', fallback=utils.sentinel):
        """Find tasks by matching the pattern with the task ID. If no matches
        are found, and fallback is not set, a ValueError will be raised. If
        fallback is set, it will be returned when no matches are found.
        Regex patterns are used by default. Set format to glob for glob
        patterns. """
        log.debug(f'Find({style}): {pattern}')

        if style == 'regex':
            regex = compile(pattern)
            matches = set([t for t in self if regex.match(t.name)])
        elif style == 'glob':
            matches = set([t for t in self if glob.fnmatch.fnmatch(t.name, pattern)])
        else:
            msg = f'Unrecognized pattern style: {style}'
            raise ValueError(msg)

        if matches:
            return matches
        elif fallback is utils.sentinel:
            err = f'No task names match value: {pattern}'
            raise ValueError(err)
        else:
            return fallback

    def graph(self):
        return WorkflowGraph(self)

    def new_task(self, *args, **kwargs):
        task = Task(*args, **kwargs)
        self.add(task)
        return task

    def pop(self, task_name):
        return self.tasks.pop(task_name)

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
        for task in self:
            task.reset()

    def resume(self):
        """Resets state for any "pending" tasks """
        log.info('Resume: Resetting state for any pending tasks...')
        for task in self:
            if task.status == 'pending':
                task.reset()

    def retry(self):
        """Resets state for any "pending" or "failed" tasks """
        log.info('Retry: Resetting state for any pending or failed tasks...')
        for task in self:
            if task.status in ('pending', 'failed', 'skipped'):
                task.reset()

    def save(self, path=None):
        save_workflow(self, path or self.path)

    def summary(self):
        return dict(Counter((t.status for t  in self)))



class WorkflowGraph:
    """Tasks are stored in separate lists in order to reduce the search time
    for the next available task. If a change to the graph occurs during
    iteration, these lists should be recalculated. """
    def __init__(self, workflow):
        log.info('Building workflow graph...')
        self.workflow = workflow
        self.G = nx.DiGraph()
        self.nodes = self.G.nodes
        self.edges = self.G.edges

        for name, task in workflow.tasks.items():
            self.G.add_node(name)

        for name, task in workflow.tasks.items():
            try:
                self._make_edges(task)
            except ValueError as e:
                err = f'While adding edges for: {task}\n{e}'
                raise ValueError(err) from None

        if not nx.is_directed_acyclic_graph(self.G):
            cycles = nx.find_cycle(self.G)
            raise ValueError(f'Not a DAG! Possible causes:{list(cycles)}')

    def __iter__(self):
        return WorkflowGraphIterator(self)

    def _add_edge(self, f, t):
        """Edges represent dependencies between tasks. Edges run FROM one node
        TO another dependent node. Nodes can have multiple edges, but not
        multiple instances of the same edge (multigraph).

            Parent ---- Is a dependency of --->  Child
         (from_node)                           (to_node)

        This means that the in-degree of a node represents the number of
        dependencies it has. A node with zero in-edges is a "root" node, or a
        task with no dependencies. """
        log.debug(f'Adding edge: {f} -> {t}')
        if f == t:
            return

        if f not in self.G:
            err = f'"{f}" is not in the workflow!'
            raise ValueError(err)

        if t not in self.G:
            err = f'"{f}" is not in the workflow!'
            raise ValueError(err)

        self.G.add_edge(f, t)


    def _make_edges(self, task):
        """Generate edges based on the floww directives of a task.

             after: task <------  target
            before: task  ------> target
             input: task <------  target

        Note: output directives do not create edges but serve as the targets of
        the input directives for other tasks
        """
        log.debug(f'Adding edges for {task}')

        for name in task.directives['after']:
            self._add_edge(name, task.name)

        for name in task.directives['before']:
            self._add_edge(task.name, name)

        for file in task.directives['input']:
            for other_task in self.workflow:
                if file in other_task.directives['output']:
                    self._add_edge(other_task.name, task.name)

        for pattern in task.directives['after-re']:
            pattern = compile(pattern)
            for other_task in self.workflow:
                if pattern.match(other_task.name):
                    self._add_edge(other_task.name, task.name)

        for pattern in task.directives['before-re']:
            pattern = compile(pattern)
            for other_task in self.workflow:
                if pattern.match(other_task.name):
                    self._add_edge(task.name, other_task.name)

        for pattern in task.directives['input-re']:
            pattern = compile(pattern)
            for other_task in self.workflow:
                for output in other_task.directives['output']:
                    if pattern.match(output):
                        self._add_edge(other_task.name, task.name)

    def ancestors(self, task):
        for anc in nx.ancestors(self.G, task.name):
            yield self.workflow[anc]

    def predecessors(self, task):
        for pre in self.G.predecessors(task.name):
            yield self.workflow[pre]

    def descendants(self, task):
        for dep in nx.descendants(self.G, task.name):
            yield self.workflow[dep]

    def successors(self, task):
        for suc in self.G.successors(task.name):
            yield self.workflow[suc]

    def is_ready(self, task):
        if task.status != 'new':
            return False

        for dep in self.predecessors(task):
            if not dep.is_complete():
                return False
        else:
            return True

    def skip_descendants(self, task, *args, **kwargs):
        for dep in self.descendants(task):
            dep.skip(dependency_failed=task.name)


class WorkflowGraphIterator:
    def __init__(self, graph):
        self.graph = graph
        self.tasks = list(self.graph.workflow)
        self.i = 0

    def __iter__(self):
        return self

    def __next__(self):
        while self.i < len(self.tasks):
            task = self.tasks[self.i]
            self.i += 1

            if self.graph.is_ready(task):
                task.pending()
                return task

        self.tasks = [t for t in self.tasks if not t.is_done()]

        if self.tasks:
            self.i = 0
            return None
        else:
            raise StopIteration


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
    wf.add(task)
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
            output = task.directives.get('output')
            inputs.append(output)

        try:
            task = jetstream.tasks.random_task(input=inputs)
            wf.add(task)
            queue.append(task)
            log.info('{} tasks added!'.format(added))
            added += 1
        except Exception as e:
            log.exception(e)

    return wf


def compile(pattern):
    return re.compile('^{}$'.format(pattern))


def load_workflow(path):
    with open(path, 'rb') as fp:
        wf = pickle.load(fp)
    wf.path = path
    return wf


def save_workflow(workflow, path):
    """Save a workflow to the path"""
    log.debug('Saving workflow: {}'.format(path))

    start = datetime.now()
    lock_path = path + '.lock'

    with open(lock_path, 'wb') as fp:
        pickle.dump(workflow, fp)

    shutil.move(lock_path, path)
    elapsed = datetime.now() - start
    log.debug('Workflow saved (after {}): {}'.format(elapsed, path))


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
    log.info(f'Mashing G:{G.path}:{len(G)} tasks with H:{H.path}:{len(H)} tasks')
    tasks = [task.copy() for task in G]
    workflow = jetstream.Workflow(tasks, props=G.props.copy())
    new = set()
    modified = set()

    for task in H:
        log.debug(f'Checking {task}')
        if task in G:
            log.debug(f'also exists in G')
            g_task = G[task.name]
            if g_task.is_failed():
                log.debug(f'status is failed, replacing in workflow..')
                workflow.pop(task.name)
                workflow.add(task)
                modified.add(task)
            elif task.identity != g_task.identity:
                log.debug(f'but identity is different, replacing in workflow..')
                workflow.pop(task.name)
                workflow.add(task)
                modified.add(task)
            else:
                log.debug(f'and same identity so skipping...')
        else:
            log.debug(f'not in G, just adding to workflow...')
            workflow.add(task)
            new.add(task)

    log.debug('Identifying tasks that need to be reset...')
    aff = new.union(modified)
    to_reset = set()
    graph = workflow.graph()

    for task in aff:
        if task.name in graph.G:
            for d in nx.descendants(graph.G, task.name):
                to_reset.add(d)
                graph.G.remove_node(d)

    for t in to_reset:
        workflow[t].reset()

    log.info(
        'Mash report:\n'
        f'New tasks: {len(new)}\n'
        f'Modified tasks: {len(modified)}\n'
        f'Total reset: {len(to_reset)}'
    )

    return workflow
