# Installation

## Install Python3 and Pip

> TGen users on Dback can load the latest version of Python with `module load python`
> and skip to the next step.

This is a Python package and requires Python3. Installation guides for Mac/Windows/Linux 
can be found here http://docs.python-guide.org/en/latest/starting/installation/
After Python3 is installed, you can install Jetstream with Pip (next step).

## Install Jetstream with Pip

Install via HTTPS (requires entering Github username and password)

```shell
pip3 install --upgrade --user git+https://github.com/tgen/jetstream.git@master
```

Install via SSH (requires SSH keys to be configured in Github profile)

```shell
pip3 install --upgrade --user git+ssh://git@github.com/tgen/jetstream.git@master
```

# Usage

Jetstream is a Python module as well as a command-line utility. 

## Command-line utility

The installation process will make two commands available:

`jetstream`

Provides access to utility/core functions for building projects and workflows.

`jetstream_pipelines`

Runs the built-in analysis pipelines.

View the help with `-h/--help` to get started. If you receive an error that a
command was not found, the Python packages bin location is probably not added 
to your `$PATH`. [Refer to this post][1] for more help.

[1]: (https://stackoverflow.com/questions/35898734/
pip-installs-packages-successfully-but-executables-not-found-from-command-line)

## Command-line Scripts

In addition to the two commands mentioned above, individual analysis modules are 
available to run as scripts and will be added to your `$PATH`. For example,
`make_exome_refkit.py -h` can be launched anywhere with no additional parameters
needed. 


## Python API

`import jetstream`

TODO: Build and link module docs, maybe readthedocs and github.io?

# Contributing

## Development install

Best practices for development are to clone the source repo, create a virtual
environment and install as a editable link with Pip. This will cause (most)
changes made to the source code to be testable without resintalling. Here is an 
example command list:

```shell
git clone https://github.com/tgen/jetstream.git
virtualenv jetstream_venv
source jetstream_venv/bin/activate
pip3 install -e jetstream/
```

## Concepts

_Project_

A directory that contains a `.jetstream` index directory. Projects can be
initialized by `jetstream project init`.

_Run_

A single instance of the jetstream application operating on a project. Each
run creates a new record in the project index directory. Each run is uniquely 
identified by its ID, a 26 character string.

_Index_

The `.jetstream` directory, and its contents, comprise a jetstream index.

_Workflow_

The set of tasks that will be executed during a run. This is modeled as a
directed acyclic graph where nodes represent tasks to complete and edges
represent dependencies between those tasks.

_Config_

Config file, run config, etc.. These are text documents describing data and
settings for the project. Project configuration files can be a variety of
text formats (csv, tsv, json, yaml) and are used for rendering workflow 
templates. They must be present in the root of the project directory, and 
will be referenced by their filename minus the extension. For example, if 
a project contained a file: `samples.csv`, it would be accessible as an array
of records via `project.samples` when rendering a workflow template. 

_Record_

[Records](https://en.wikipedia.org/wiki/Record_(computer_science)) are 
generated for any data stored in tables (csv or tsv). Essentially, each row 
in the table becomes an object with key-value properties determined by the 
header of the table. In Python, the table becomes a list of dictionaries
available in the project object. This strategy allows project data to be 
created/stored in a ton of different formats, while needing only a single 
interface to access that data (see the Project() class).






