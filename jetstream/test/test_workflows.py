import jetstream
import logging



def test_workflow():
    logging.root.setLevel(logging.DEBUG)
    cmd = 'echo hello world'

    wf = jetstream.Workflow()
    wf.add_node('first_task', cmd=cmd)
    wf.add_node('second_task', cmd=cmd)
    wf.add_node('third_task', cmd=cmd)
    wf.add_node('start', cmd=cmd)
    wf.add_node('end', cmd=cmd)

    wf.add_dependency('start', before='.*_task')
    wf.add_dependency('end', after='.*_task')

    return wf
