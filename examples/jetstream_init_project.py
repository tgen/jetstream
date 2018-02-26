#!/usr/bin/env python3
import os
import sys
import jetstream

# Load the config file given as first argument
config = jetstream.config.legacy.load(sys.argv[1])

# Create a project dir
proj_dir = config['meta']['project']
os.mkdir(proj_dir)
os.mkdir(os.path.join(proj_dir, '.jetstream'))

# Save the project data inside the project dir
config_path = os.path.join(proj_dir, 'project.json')
with open(config_path, 'w') as fp:
    fp.write(jetstream.config.serialize(config, format='json'))

print('Initialized project {}'.format(config_path))
