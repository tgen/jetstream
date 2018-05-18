import jetstream
import logging


def test_workflow():
    logging.root.setLevel(logging.DEBUG)
    cmd = 'echo hello world'

    wf = jetstream.Workflow()
    wf.add_task('first_task', cmd=cmd)
    wf.add_task('second_task', cmd=cmd)
    wf.add_task('third_task', cmd=cmd)
    wf.add_task('start', cmd=cmd)
    wf.add_task('end', cmd=cmd)

    wf.add_dependency('start', before='.*_task')
    wf.add_dependency('end', after='.*_task')

    return wf
