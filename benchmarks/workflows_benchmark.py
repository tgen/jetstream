#!/usr/bin/env python3
import sys
import jetstream

print('Starting!')

wf = jetstream.Workflow()

for i in range(int(sys.argv[1])):
    wf.add_task(str(i))

print('Workflow built: {} tasks'.format(len(wf)))

for task in wf.tasks():
    wf.complete(task)

print('Workflow complete!')

