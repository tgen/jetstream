import json
import logging
import shutil
from datetime import datetime

import networkx as nx
from networkx.drawing.nx_pydot import to_pydot
from networkx.readwrite import json_graph

from jetstream import utils, exc


log = logging.getLogger(__name__)


class Workflow:
    """ Workflows are a network graph representating a series of steps that need
    to be executed on a project. Each node of the graph represents a plugin
    component loaded with jetstream.plugins.get_plugin(). Workflows are a type
    of iterable that implement __next__() and __send__() methods. Calling next()
    on a workflow will yield a plugin that is ready for to be executed, or None
    if there are none available. StopIteration will be raised when all plugins
    are completed.

    on_update takes a function that will be called every time a change is made
    to the workflow.

    Workflows can be loaded from existing graphs with the graph argument or
    from_pydot() method.

    """

    def __init__(self, *, data=None, **attr):
        self.last_update = dict()
        self.graph = data or nx.DiGraph(**attr)

        if not nx.is_directed_acyclic_graph(self.graph):
            raise exc.NotDagError

    # Utility methods
    def __str__(self):
        """ Gives better results when using print() """
        return utils.yaml_dumps(self.serialize())

    def save(self, path):
        lock_path = path + '.lock'
        obj = {'workflow': self.serialize()}
        utils.yaml_dump(obj, lock_path)
        shutil.move(lock_path, path)

    def serialize(self, serializer=None):
        """Serialize the workflow for dumping to json/yaml etc..."""
        if serializer is not None:
            return serializer(self.graph)
        else:
            return self.to_node_link_data()

    def to_dot(self):
        """ Returns a pydot representation of the graph """
        # Note Had to dig into the networkx source to find this, but it works..
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
        # TODO we may want this to reset the status for pending tasks?
        return self

    def __next__(self):
        """ returns the next component available for launch and marks
        the node status as "pending" """
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

        self.update(node_id, datetime_end=datetime.now())

        if return_code != 0:
            self.fail(node_id, return_code=return_code, logs=logs)
        else:
            self.complete(node_id, return_code=return_code, logs=logs)

    def complete(self, node_id, return_code=0, logs=''):
        """ Complete a node_id """
        log.critical('Node complete! {}'.format(node_id))
        self.update(
            node_id,
            status='complete',
            return_code=return_code,
            logs=logs
        )

    def fail(self, node_id, return_code=1, logs=''):
        """ Fail a node_id, this also fails any nodes dependent on node_id"""
        log.critical('Node failed! {}'.format(node_id))

        self.update(
            node_id,
            status='failed',
            return_code=return_code,
            logs=logs
        )

        for d in self.graph.predecessors(node_id):
            self.fail(d, logs='Failed due to dependency {}'.format(node_id))

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
            raise exc.NotDagError

    def add_node(self, id, cmd, **kwargs):
        log.debug('')
        node_data = {
            'id': id,
            'cmd': cmd,
        }

        node_data.update(**kwargs)
        self._add_node(id, node_data)

        return id

    def add_dependency(self, id, before=None, after=None, **kwargs):
        if before:
            if isinstance(before, str):
                before = (before,)

            for child_id in before:
                if not child_id in self.graph:
                    raise ValueError('{} not in graph'.format(child_id))
                self._add_edge(from_node=child_id, to_node=id)

        if after:
            if isinstance(after, str):
                after = (after,)

            for parent_id in after:
                if not parent_id in self.graph:
                    raise ValueError('{} not in graph'.format(parent_id))
                self._add_edge(from_node=id, to_node=parent_id)

    def compose(self, wf):
        res = nx.algorithms.binary.compose(self.graph, wf.graph)
        self.graph = res

#
# class Result(object):
#     """A common structure for the launchers to return to the workflow"""
#
#     def __init__(self, plugin, logs, return_code, error=None):
#         self.fingerprint = utils.fingerprint()
#         self.plugin = plugin
#         self.logs = str(logs)
#         self.return_code = int(return_code)
#         self.error = str(error)
#
#     def to_json(self):
#         return json.dumps(self.serialize())
#
#     def serialize(self):
#         """Returns a dictionary ready for json serialization. """
#         return {
#             'plugin': self.plugin,
#             'fingerprint': self.fingerprint,
#             'logs': str(self.logs),
#             'return_code': int(self.return_code),
#             'error': str(self.error),
#
#         }


def from_node_link_data(data):
    graph = json_graph.node_link_graph(data)
    return Workflow(data=graph)


def to_node_link_data(wf):
    return json_graph.node_link_data(wf.graph)
