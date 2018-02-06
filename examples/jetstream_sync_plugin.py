#!/usr/bin/env python3
import os
import sys
import jetstream

plugin = sys.argv[1]

if plugin.startswith(('https', 'git@github.com')):
  pass
else:
  plugin = os.path.realpath(plugin)
  
jetstream.plugins._clone(plugin)
