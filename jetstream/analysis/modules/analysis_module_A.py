"""This is the placeholder for an analysis module. These modules can be any application/script that is platform-compatible
with the API-Command stack (RHEL in our case). They communicate logs back to the command module through stderr/stdout/
The command module will decide whether or not to continue with the workflow based on the exit code, so make it meaningful.

These applications will always be started in a project directory that looks like this:
    myproject-1234/
    ├── data
    │   └── project.json
    └── tmp

They have read-only access in data/, and read-write access in tmp/. A project manifest will always be available in data/ that 
stores project-wide runtime variables. These variables can be used by the analysis modules for configuration. 
"""
import time
import os

print('Analysis Module A is running on: ' + os.getcwd())
time.sleep(300)
print("Analysis Module A is done!")
