#!/usr/bin/env python3
import sys
from datetime import datetime
import jetstream

print('Starting!')
start = datetime.now()

wf = jetstream.Workflow()

for i in range(int(sys.argv[1])):
    wf.new_task(str(i))

print('Workflow built: {}'.format(wf))

i = iter(wf)

try:
    while 1:
        t, d = next(i)
        i.send(t, 0)
except StopIteration:
    print('StopIteration')
finally:
    print(datetime.now() - start)

print(wf)
print('Workflow complete!')

