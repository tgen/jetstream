from jetstream.core import workflows


def workflow1():
    wf = workflows.Workflow()
    wf.add_node('A', ['env'])
    wf.add_node('B', ['env'])
    wf.add_node('C', ['hostname'])
    return wf
