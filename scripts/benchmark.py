#!/usr/bin/env python3

"""A very simple benchmarking tool."""

import sys
import timeit
import subprocess
import statistics

def squelch_call():
    subprocess.call(
        sys.argv[1:],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL)

t = timeit.Timer(stmt=squelch_call)

res = []
trial = 0
while len(res) < 1000:
    trial += 1
    if trial % 100 == 0:
        print('Trial:', trial, file=sys.stderr)
 
    try:
        res.append(t.timeit(1))
    except KeyboardInterrupt:
        break

print('\nn Trials:', len(res))
print('min:', min(res))
print('max:', max(res))
print('mean:', statistics.mean(res))
print('median:', statistics.median(res))
print('stdev:', statistics.stdev(res))
