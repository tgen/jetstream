
# Installation

## Make sure you have Python3 installed

If you already have Python3 and Pip installed, skip to the next step.

This requires Python3. For a lot of reasons, [Homebrew](https://brew.sh) is 
the best way to install Python3. Follow this guide:

http://docs.python-guide.org/en/latest/starting/install3/osx/

## Single Command Install with Pip

Install via HTTPS (requires entering Github username and password)

```shell
pip3 install --upgrade --user git+https://github.com/tgen/jetstream.git@master
```

Install via SSH (requires SSH keys to be configured in Github profile)

```shell
pip3 install --upgrade --user git+ssh://git@github.com/tgen/jetstream.git@master
```

# Development notes

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

The set of tasks that will be executed during a run. This is organized as a
directed acyclic graph to respect plugin dependencies.

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
header of the table. This strategy allows project to be stored in a ton of 
different formats, while needing only a single interface to access that 
data (see the Project() class).







