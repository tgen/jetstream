import jetstream

wf = jetstream.Workflow()
wf.add_task('first_task', cmd='echo hello world')
wf.add_task('second_task', cmd='env')
wf.add_task('third_task', cmd='hostname')
wf.add_task('start', cmd=None)
wf.add_task('end', cmd=None)

wf.add_dependency('start', before='.*_task')
wf.add_dependency('end', after='.*_task')
