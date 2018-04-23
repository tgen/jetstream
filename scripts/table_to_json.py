#!/usr/bin/env python3
import sys
import json
import jetstream

print(json.dumps(jetstream.utils.table_to_records(sys.argv[1])))
