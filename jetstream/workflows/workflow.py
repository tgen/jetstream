import json
import re
import logging
import shutil
from datetime import datetime
import shlex
import networkx as nx
from networkx.drawing.nx_pydot import to_pydot
from networkx.readwrite import json_graph

from jetstream import utils

log = logging.getLogger(__name__)


class NotDagError(Exception):
    """ Raised when an action would result in a network that is not a DAG """
    pass


class Workflow:
    """ Workflows are a network graph representing a series of steps that need
    to be executed on a project. Each node of the graph represents a task
    to be completed, edges represent dependencies between the tasks. Workflows
    are iterables that implement __next__() and __send__() methods. Calling
    next() on a workflow will return a task (node_id, node_data) that is ready
    for to be executed, or None if there are none currently available.
    StopIteration will be raised when all tasks are completed. Task status
    should be updated via Workflow.__send__(), Workflow.fail_node(), or
    Workflow.pass_node().

    Workflows can be loaded from existing graphs with the data argument or
    from_node_link_data() method.

    """

    def __init__(self, *, data=None, **attr):
        self.last_update = dict()
        self.graph = data or nx.DiGraph(**attr)

        if not nx.is_directed_acyclic_graph(self.graph):
            raise NotDagError

    def __str__(self):
        """ Gives better results when using print() """
        return utils.yaml_dumps(self.serialize())

    def serialize(self, serializer=None):
        """Serialize the workflow for dumping to json/yaml etc..."""
        if serializer is not None:
            return serializer(self.graph)
        else:
            return self.to_node_link_data()

    def to_dot(self):
        """ Returns a pydot representation of the graph. Note: Had to dig into
        the networkx source to find this, but it works.."""
        graph = self.graph.copy()
        return to_pydot(graph).__str__()

    def to_yaml(self):
        return utils.yaml_dumps(obj=self.serialize())

    def to_json(self):
        return json.dumps(self.serialize(), indent=4)

    @staticmethod
    def from_node_link_data(*args, **kwargs):
        return from_node_link_data(*args, **kwargs)

    def to_node_link_data(self, wf=None):
        if wf is None:
            return to_node_link_data(self)
        else:
            return to_node_link_data(wf)

    # Methods for iterating over a workflow
    def __iter__(self):
        """ Workflows are essentially fancy iterators that allow feedback
        through the __send__ method. See generators """
        return self

    def __next__(self):
        """ Returns the next task available for launch and marks the node status
        as "pending" """
        pending = False
        for node_id, node_data in self.graph.nodes(data=True):
            if node_data['status'] == 'pending':
                pending = True
            if self.node_ready(node_id):
                log.debug('Releasing node for execution: {}'.format(node_id))

                # Mark the node as pending
                self.update(
                    node_id,
                    status='pending',
                    datetime_start=str(datetime.now())
                )

                return node_id, node_data
        else:
            if pending:
                return None
            else:
                log.debug('Request for next task but all complete')
                raise StopIteration

    def __send__(self, node_id, return_code, logs):
        """ Returns results to the workflow """
        log.debug('Received results for {}'.format(node_id))

        self.update(node_id, datetime_end=str(datetime.now()))

        if return_code != 0:
            self.fail_node(node_id, return_code=return_code, logs=logs)
        else:
            self.pass_node(node_id, return_code=return_code, logs=logs)

    def pass_node(self, node_id, return_code=0, logs=''):
        """ Complete a node_id """
        log.critical('Node complete! {}'.format(node_id))
        self.update(
            node_id,
            status='complete',
            return_code=return_code,
            logs=logs
        )

    def fail_node(self, node_id, return_code=1, logs=''):
        """ Fail a node_id, this also fails any nodes dependent on node_id"""
        log.critical('Node failed! {}'.format(node_id))

        self.update(
            node_id,
            status='failed',
            return_code=return_code,
            logs=logs
        )

        for d in self.graph.predecessors(node_id):
            self.fail_node(d, logs='Failed due to dependency {}'.format(node_id))

    def reset_node(self, node_id):
        log.critical('Node reset! {}'.format(node_id))

        self.update(
            node_id,
            status='new',
            return_code=None,
            logs=None,
            datetime_start=None,
            datetime_end=None
        )

    # Methods for reading from a workflow
    def nodes(self, *args, **kwargs):
        return self.graph.nodes(*args, **kwargs)

    def status(self, node_id=None):
        """ Returns the status of a node or all nodes if node is None """
        if node_id is not None:
            return self.graph.nodes[node_id]['status']
        else:
            return [d['status'] for n, d in self.graph.nodes(data=True)]

    def get_node(self, node_id, data=False):
        """ Returns None if the node is not in the graph """
        g = self.graph
        if node_id in g:
            if data:
                return g.nodes()[node_id]
            else:
                return node_id
        else:
            return None

    def node_ready(self, node_id):
        """ Returns True if the given node is ready for execution """
        g = self.graph
        node_data = g.nodes()[node_id]

        if node_data['status'] != 'new':
            return False

        for dependency in self.graph.successors(node_id):
            dependency_status = g.nodes()[dependency]['status']
            if dependency_status != 'complete':
                return False
        else:
            return True

    def root_nodes(self):
        """ Returns the set of root nodes in the graph """
        graph = self.graph

        res = set()
        for node in graph.nodes():
            # Note: inspectors warn about calls to graph.out_degree (probably)
            # because it's an overloaded method. But this line seems to work
            # as intended. The goal is to find nodes in the graph with an
            # out_degree of 0.
            if graph.out_degree(node) == 0:
                res.add(node)
        return res

    def update(self, node_id, **kwargs):
        """ Change the status of a node """
        log.debug('Update node {}: {}'.format(node_id, kwargs))
        self.last_update = utils.fingerprint()
        self.graph.nodes[node_id].update(**kwargs)

    def _add_node(self, node_id, data):
        """ Adding a node requires a mapping that includes "name" key. A
        RuntimeError will be raised if the name already exists in the graph"""

        if self.get_node(node_id):
            raise RuntimeError('Duplicate node id: {} in\n{}'.format(
                id, self))

        # All nodes get a status attribute that the workflow uses to identify
        # nodes that are ready to be executed.
        data['status'] = 'new'
        return self.graph.add_node(node_id, **data)

    def _add_edge(self, from_node, to_node):
        """ Edges represent dependencies between components. Edges run
        FROM one node TO another node that it depends upon. Nodes can have
        multiple edges, but not multiple instances of the same edge.

            Child ----- Depends Upon -----> Parent
         (from_node)                       (to_node)

        This means that the out-degree of a node represents the number
        of dependencies it has. A node with zero out-edges is a "root"
        node, or a component with no dependencies.

        The methods ".add_component_before()" and ".add_component_after()" are
        provided for adding components to a workflow, and should be preferred
        over adding edges directly to the workflow.
        """
        g = self.graph

        g.add_edge(from_node, to_node)

        if not nx.is_directed_acyclic_graph(g):
            g.remove_edge(from_node, to_node)
            raise NotDagError

    def add_node(self, node_id, **kwargs):
        if 'cmd' not in kwargs or not kwargs['cmd']:
            msg = "Node: '{}' is missing a required value for 'cmd'"
            raise ValueError(msg.format(node_id))

        if isinstance(kwargs['cmd'], str):
            kwargs['cmd'] = shlex.split(kwargs['cmd'])

        self._add_node(node_id, kwargs)
        return node_id

    def add_dependency(self, node_id, before=None, after=None):
        if not node_id in self.graph:
            raise ValueError('{} not in graph'.format(node_id))

        stack = []
        try:
            if before:
                if isinstance(before, str):
                    before = (before,)

                before_pats = [re.compile(b) for b in before]

                for pat in before_pats:
                    matches = filter(pat.match, self.nodes())

                    if not matches:
                        raise ValueError('No matching nodes for: {}'.format(
                            pat.pattern))

                    for child_id in matches:
                        self._add_edge(from_node=child_id, to_node=node_id)
                        stack.append((child_id, node_id))

            if after:
                if isinstance(after, str):
                    after = (after,)

                after_pats = [re.compile(a) for a in after]

                for pat in after_pats:
                    matches = filter(pat.match, self.nodes())

                    if not matches:
                        raise ValueError('No matching nodes for {}'.format(
                            pat.pattern))

                    for parent_id in matches:
                        self._add_edge(from_node=node_id, to_node=parent_id)
                        stack.append((node_id, parent_id))

        except Exception as e:
            for u, v in stack:
                log.debug('Rolling back: {} -> {}'.format(u, v))
                self.graph.remove_edge(u, v)
            raise e from None

    def compose(self, wf):
        res = nx.algorithms.binary.compose(self.graph, wf.graph)
        self.graph = res

    def resume(self):
        """ Returns all pending nodes to an incomplete state. """
        for node_id, node_data in self.nodes(data=True):
            if node_data['status'] == 'pending':
                self.reset_node(node_id)

    def reset(self):
        """ Returns all nodes to a new state. """
        for node_id, node_data in self.nodes(data=True):
            self.reset_node(node_id)


def from_node_link_data(data):
    graph = json_graph.node_link_graph(data)
    return Workflow(data=graph)


def to_node_link_data(wf):
    return json_graph.node_link_data(wf.graph)


def load(path, loader=utils.yaml_load):
    data = loader(path)
    return from_node_link_data(data['workflow'])


def save(wf, path):
    lock_path = path + '.lock'
    obj = {'workflow': wf.serialize()}
    utils.yaml_dump(obj, lock_path)
    shutil.move(lock_path, path)
