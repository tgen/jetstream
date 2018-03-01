import json
import logging
import networkx as nx
from datetime import datetime
from uuid import uuid4 as uuid
from networkx.drawing.nx_pydot import to_pydot
from networkx.readwrite import json_graph

from jetstream import plugins, utils

log = logging.getLogger(__name__)


class NotDagError(Exception):
    """ Raised when an action would result in a network that is not a DAG """
    pass


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
    def __init__(self, *, on_update=None, graph=None):
        self.on_update = on_update or log.debug
        self.graph = graph or nx.DiGraph()

        if not nx.is_directed_acyclic_graph(self.graph):
            raise NotDagError

    # Utility methods
    def __str__(self):
        """ Gives better results when using print() """
        return utils.yaml_dump(self.serialize())

    def serialize(self, serializer=None):
        """Serialize the workflow for dumping to json/yaml etc..."""
        if serializer is not None:
            return serializer(self.graph)
        else:
            data = nx.readwrite.json_graph.node_link_data(self.graph)
            return data

    def to_dot(self):
        """ Returns a pydot representation of the graph """
        # Note Had to dig into the networkx source to find this, but it works..
        graph = self.graph.copy()
        return to_pydot(graph).__str__()

    def to_yaml(self, *args, **kwargs):
        return utils.struct(action='dumps', format='yaml', obj=self.serialize())

    def to_json(self, *args, **kwargs):
        return utils.struct(action='dumps', format='json', obj=self.serialize())

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
                log.debug('Node ready for execution {}'.format(node_id))
                self.update(
                    node_id,
                    status='pending',
                    datetime_start=str(datetime.now())
                )
                return node_id, node_data
        else:
            if pending:
                log.debug('Request for next task but None available')
                return None
            else:
                log.debug('Request for next task but all complete')
                raise StopIteration

    def __send__(self, node_id, result):
        """ Returns results to the workflow """
        log.critical('Received results for {}'.format(node_id))

        if result.return_code != 0:
            # TODO handle results better
            # this needs to recognized failures and set node status to new
            # but we might also want to only allow a limited number of retrys
            # per node or globally
            log.critical('Node returned non-zero, here we would decide'
                         'how many times to try and repeat that node, or'
                         'halt and inform the user. ')

        self.update(
            node_id,
            status='complete',
            datetime_end=datetime.now(),
            **result.serialize()
        )

    # Methods for reading from a workflow
    def nodes(self, *args, **kwargs):
        return self.graph.nodes(*args, **kwargs)

    def components(self, *args, **kwargs):
        """ alias for Workflow.nodes() """
        return self.nodes(*args, **kwargs)

    def status(self, node=None):
        """ Returns the status of a node or all nodes if node is None """
        if node is not None:
            return self.graph.nodes[node]['status']
        else:
            return [d['status'] for n, d in self.graph.nodes(data=True)]

    def get_node(self, name, data=False):
        """ Returns None if the node is not in the graph """
        g = self.graph
        if name in g:
            if data:
                return g.nodes()[name]
            else:
                return name
        else:
            return None

    def node_ready(self, node):
        """ Returns True if the given node is ready for execution """
        g = self.graph
        node_data = g.nodes()[node]

        if node_data['status'] != 'new':
            return False

        for dependency in self.graph.successors(node):
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
            # Note: inspectors warn about calls to graph.out_degree, probably
            # because it's an overloaded method. But this line seems to work
            # as intended. We're trying to find nodes in the graph with an
            # out_degree of 0.
            if graph.out_degree(node) == 0:
                res.add(node)
        return res

    def update(self, node, **kwargs):
        """ Change the status of a node """
        self.last_update = utils.fingerprint()
        self.graph.nodes[node].update(**kwargs)

    def _add_node(self, mapping):
        """ Adding a node requires a mapping that includes "id" key. The id will
        be used to reference the node in the graph, the rest of the mapping will
        be added as data """
        node_id = str(uuid())
        self.graph.add_node(node_id, **mapping)
        return self.get_node(node_id)

    def _add_edge(self, from_node, to_node):
        """ Edges represent dependencies between components. Edges run
        FROM one node TO another node that the component depends upon. Nodes
        can have multiple edges.

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

    def add_component_before(self, plugin_id, *before):
        """ Add a component and specify that it should run before some other
        component(s) """
        for child_id in before:
            if not child_id in self.nodes():
                raise ValueError('Node: {} not in workflow'.format(child_id))
        else:
            parent_id = self.add_component(plugin_id)
            for node_id in before:
                child_id = self.get_node(node_id)
                self._add_edge(from_node=child_id, to_node=parent_id)
            return parent_id

    def add_component_after(self, plugin_id, *after):
        """ Add a component and specify that it should run after some other
        component(s) """
        for parent_id in after:
            if not parent_id in self.nodes():
                raise ValueError('Node: {} not in workflow'.format(parent_id))
        else:
            child_id = self.add_component(plugin_id)
            for parent_id in after:
                self._add_edge(from_node=child_id, to_node=parent_id)
            return child_id

    def add_component(self, plugin_id):
        plugin = plugins.get_plugin(plugin_id)
        node_data = dict(status='new', **plugin)
        return self._add_node(node_data)

    def add_component_from_file(self, path):
        plugin = plugins.load(path)
        node_data = dict(status='new', **plugin)
        return self._add_node(node_data)


class Result(object):
    """A common structure for the launchers to return to the workflow"""
    def __init__(self, plugin_id, log, return_code, error=None):
        self.plugin_id = str(plugin_id)
        self.log = str(log)
        self.return_code = int(return_code)
        self.error = str(error)

    def to_json(self):
        return json.dumps(self.serialize())

    def serialize(self):
        """Returns a dictionary ready for json serialization. """
        return {
            'plugin_id': str(self.plugin_id),
            'log': str(self.log),
            'return_code': int(self.return_code),
            'error': str(self.error)
        }


def from_node_link_data(data, *args, **kwargs):
    graph = json_graph.node_link_graph(data)
    return Workflow(graph=graph, *args, **kwargs)


def to_node_link_data(wf):
    return json_graph.node_link_data(wf.graph)

