from .workflow import Workflow, serialize_json, serialize_pydot

def load_pydot(path):
    from pydot import graph_from_dot_file
    graph = graph_from_dot_file(path)[0]
    wf = Workflow.from_pydot(graph)
    return wf