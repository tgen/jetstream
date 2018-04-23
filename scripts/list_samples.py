#!/usr/bin/env python3
""" Prints JSON records for all samples in a project. """
import jetstream
import json

p = jetstream.Project()

for sample in p.list_samples():
    print(json.dumps(sample))
