import json
import logging
from uuid import uuid4 as uuid
import networkx as nx
import pydot
from networkx.drawing.nx_pydot import to_pydot, from_pydot

log = logging.getLogger(__name__)


class NotDagError(Exception):
    """ Raised when an action would result in a network that is not a DAG """
    pass


class Workflow:
    def __init__(self, id=None, notify=None, graph=None):
        self.id = id or str(uuid())
        self.notify = notify or log.debug
        self.graph = graph or nx.DiGraph(name=self.id, )

    # Utility methods
    def __str__(self):
        """ Gives better results when using print() """
        data = self.serialize_pydot()
        return json.dumps(data, indent=4)

    def _clean_graph(self):
        graph = self.graph.copy()
        for node in graph:
            print(node['status'])

    def serialize_json(self):
        """ Returns a json representation of the graph """
        data = nx.json_graph.cytoscape_data(self.graph)
        return data

    def serialize_pydot(self):
        # Note Had to dig into the networkx source to find this, but it works..
        graph = self.graph.copy()
        return to_pydot(self.graph).__str__()

    @staticmethod
    def from_pydot(p, *args, **kwargs):
        graph = from_pydot(p)
        return Workflow(graph=graph, *args, **kwargs)

    # Methods for iterating over a workflow
    def __iter__(self):
        """ Workflows are essentially fancy iterators that allow feedback
        through the __send__ method. See generators """
        # we may want this to reset all the tasks?
        return self

    def __send__(self, result):
        """ Returns results to the workflow """
        node_name, record = result
        log.critical('Workflow recieved results for {}\n{}'.format(node_name, record))
        # TODO handle results better
        self.update(node_name, 'complete')

    def __next__(self):
        """ returns the next module available for launch and marks its status as "pending" """
        pending = False
        for node, node_data in self.graph.nodes(data=True):
            if node_data['status'] == 'pending':
                pending = True
            if self.node_ready(node):
                self.update(node, 'pending')
                return (node, node_data)
        else:
            if pending:
                return None
            else:
                raise StopIteration

    # Access/update modules in the workflow
    def nodes(self, *args, **kwargs):
        return self.graph.nodes(*args, **kwargs)

    def modules(self, *args, **kwargs):
        """ same as Workflow.nodes() """
        return self.nodes(*args, **kwargs)

    def update(self, node, status):
        """ Change the status of a node """
        self.graph.nodes[node]['status'] = status

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
            if graph.out_degree(node) == 0:
                res.add(node)
        return res

    def _auto_notify(callback):
        """ This decorator sends a tuple of useful data to Workflow.notify()
        prior to the function being called """
        def fn(self, *args, **kwargs):
            if self.notify is not None:
                self.notify((self.id, callback.__name__, args, kwargs))
            return callback(self, *args, **kwargs)
        return fn

    @_auto_notify
    def _add_node(self, name, **kwargs):
        self.graph.add_node(name, **kwargs)
        return self.get_node(name)

    @_auto_notify
    def _add_edge(self, from_node, to_node):
        """ Edges represent dependencies between modules. Edges run
        FROM one node TO another node that the module depends upon. Nodes
        can have multiple edges.

            Child ----- Depends Upon -----> Parent
         (from_node)                       (to_node)

        This means that the out-degree of a node represents the number
        of dependencies it has. A node with zero out-edges is a "root"
        node, or a module with no dependencies.

        The methods ".add_module_before()" and ".add_module_after()" are
        provided for adding modules to a workflow, and should be preferred 
        over adding edges directly to the workflow.
        """ 
        g = self.graph

        g.add_edge(from_node, to_node)

        if not nx.is_directed_acyclic_graph(g):
            g.remove_edge(from_node, to_node)
            raise NotDagError

    @_auto_notify
    def add_module_before(self, before, *names):
        """ Add a module and specify that it should run before some other
        module(s) """

        child = self.get_node(before) or self.add_module(before)
        for n in names:
            parent = self.get_node(n) or self.add_module(n)
            self._add_edge(from_node=child, to_node=parent)

    @_auto_notify
    def add_module_after(self, after, *names):
        """ Add a module and specify that it should run after some other
        module(s) """
        parent = self.get_node(after) or self.add_module(after)

        for n in names:
            child = self.get_node(n) or self.add_module(n)
            self._add_edge(from_node=child, to_node=parent)

    def add_module(self, module):
        self._add_node(module, status='new')
        return module

