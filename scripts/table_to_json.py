#!/usr/bin/env python3
import sys
import json
import jetstream

print(json.dumps(jetstream.utils.load_table(sys.argv[1])))
