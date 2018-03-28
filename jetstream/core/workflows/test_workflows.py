from jetstream.core.workflows.workflow import Workflow


def workflow1():
    wf = workflows.Workflow()
    wf.add_node('A', ['env'])
    wf.add_node('B', ['env'])
    wf.add_node('C', ['hostname'])
    return wf
