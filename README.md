# Installation

## Install Python3 and Pip

> TGen users on Dback can load the latest version of Python with `module load python` and skip to the next step.

This is a Python package and requires Python3. Installation guides for Mac/Windows/Linux are available from the [Hitchiker's Guide to Python][install_help] After Python3 is installed, you can install Jetstream with Pip (next step).


## Install Jetstream with Pip

Install via HTTPS (requires entering Github username and password)

```shell
pip3 install --upgrade --user git+https://github.com/tgen/jetstream.git@master
```

Install via SSH (requires SSH keys to be configured in Github profile)

```shell
pip3 install --upgrade --user git+ssh://git@github.com/tgen/jetstream.git@master
```

# Intro

## Workflow directives

All directives are optional. It's possible (but maybe useless) to make tasks
without a `cmd`, for example. Some directives will be used differently 
depending on the backend you're running with: Local vs Slurm. Task and
directive order is irrelevant.

```yaml
- name: Name for linking flow to this task
  after: 
    - This task will run after tasks matching
    - each given value, Supports Python regex 
    - patterns
  before: 
    - This task will run before tasks matching
    - each given value. Supports Python regex 
    - patterns
  input:
    - This task will run after tasks with
    - output directives matching each
    - given value. Supports Python regex
    - patterns
  output:
    - This task will satisfy a matching
    - input value requirement.
  cmd: Shell command executed with /bin/bash
  stdin: Data will be connected to stdin of cmd
  stdout: cmd stdout will be connected to value
  stderr: cmd stderr will be connected to value
  cpus: LocalBackend - Will reserve local cpus when 
    launching cmd
    (SlurmBackend) "-c" when requesting job
    allocation
  mem: SlurmBackend - "--mem" when requesting
    job allocation
  # Other optional directives
  methods: Describe what this task does in plain
    language, later this can be used to 
    generate methods for a project.
  description: Describe this task for potential 
    users
```

# Docs

Full documentation can be found here http://dback-login3.tgen.org:8082 , they're still largely a work-in-progress:


[install_help]: http://docs.python-guide.org/en/latest/starting/installation/
[path_help]: https://stackoverflow.com/questions/35898734/pip-installs-packages-successfully-but-executables-not-found-from-command-line
