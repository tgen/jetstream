# Jetstream

Jetstream is a pipeline development framework written as a pure Python package 
and command line utility. It supports complex workflows modeled as directed-
acyclic graphs (DAGs), and execution on batch schedulers (Slurm). This document 
will introduce the core concepts of framework, and give some real examples.

## Install

> TGen users on Dback can load the latest version of Python with `module load python`.

The easiest way to install Jetstream is with Pip:

```shell
pip3 install --upgrade --user git+https://github.com/tgen/jetstream.git@master
```

More install help can be found in the [installation](#installation) section below

## Command-line

After installing with Pip, the command line utility can be launched with
`jetstream`, and help can be accessed with the `-h/--help` options. If 
the command is not found, see the detailed [installation](#installation) 
help below.

```shell
$ jetstream -h
usage: jetstream [-v] [--log-debug] [--log-verbose] [--log-format LOG_FORMAT]
                 [--log-filename LOG_FILENAME] [--log-filemode LOG_FILEMODE]
                 [--log-level LOG_LEVEL]

Available commands are:

    build
        Build workflow files.

    draw
        Draw the network graph for a workflow.

    init
        Create a new project or reinitialize an existing project.

    project
        Interact with jetstream projects. View tasks, run history, or
        project data. This option requires a subcommand.

    run
        Run Jetstream from a template, module, or workflow.

optional arguments:
  -v, --version         show program's version number and exit
  --log-debug           Alias for debug log settings
  --log-verbose         Alias for lowest-level log settings
  --log-format LOG_FORMAT
  --log-filename LOG_FILENAME
  --log-filemode LOG_FILEMODE
  --log-level LOG_LEVEL
```
To run a workflow template with a command line argument use the commands shown in the following example.
 ```
 jetstream run commandLine.jst --str:input inputFilename.txt
```
The command line argument `--str:input inputFilename.txt` follows the format `--<type>:key value`.
The `input` variable name corresponding to `key` is referenced in the `commandLine.jst` template 
within double-nested braces `{{input}}`, as described in the [***Variables in templates***](#variables-in-templates) section.

## Python package

```python
>>> import jetstream
>>> wf = jetstream.Workflow()
>>> wf.new_task(name='task1')
<Task task1>
>>> wf.save(path='mywf.yaml')
```

Full package documentation can be found here:
http://dback-login3.tgen.org:8082

## Contributing

All contributions are welcome.
Comments, suggestions, pull-requests, and issues can be submitted on [https://github.com/tgen/jetstream/issues]


# Introduction

Pipelines are a set of computational tasks that need to be repeated for multiple sets 
of input data. Jetstream models pipelines as _Tasks_ and _Workflows_:

- **Tasks** are the fundamental unit of a pipeline. They describe what to do: usually
  a shell command to execute. And when to do it: with directives like `before`, `after`,
  or `input` that establish dependencies between tasks.

  Tasks can include variables that will be filled-in by configuration data. This allows 
  a single set of tasks to be reused for multiple datasets (ie a pipeline). For each 
  run, the tasks will be combined with new configuration data to produce the finalized 
  runner instructions. A set of tasks that have been finalized with config data is 
  called a Workflow. 

- **Workflows** are the internal data structure that models a group of tasks as a DAG.
  They can be thought of as _specific_ instances of a _generalized_ set of tasks. The runner 
  uses the workflow object to coordinate the executions of the tasks, and store the progress
  of the run. Users rarely need to interact with workflows directly, but they are a critical
  step in the operation of pipelines.

There are two ways that pipelines can be built with Jetstream:

- [**Templates**](#building-pipelines-from-templates): Define the tasks in a text document
  that can be processed and run with the `jetstream` command-line utility.

- [**Python Module**s](#building-pipelines-from-python-modules): Using the functions and
  classes in the `jetstream` package to manually construct a workflow object.


# Building Pipelines from Templates

**Templates** are hybrid YAML-Jinja2 text files (.jst extension by convention) that describe 
a skeleton of the tasks in a pipeline. Once the template is combined with config data, the tasks
will be a complete description of what needs to happen for that single run of the pipeline. In
YAML terms, they are a sequence of mappings, where each mapping element represents single task, 
and each item in the mapping is a task directive. But basically, each task starts with `-` and then 
task directives are given as `directive: value` lines. YAML is a whitespace oriented file format, 
so pay attention to the spaces, and use a text editor with syntax highlighting for YAML.

Here is an example template that describes two tasks to complete.

```yaml

- name: task1
  cmd: hostname

- name: task2
  cmd: date

```

Each task starts on a line with `-`, this signals the begining of an item in the sequence. Next,
there are directives defined. In this case the tasks are relatively simple: they have a `name`
which allows us to link dependencies, and a `cmd` that will be executed.

After the template is created, it can be ran with the workflow runner. Save the
template above in a text file called 'example.jst'. Then, run it with the command
`jetstream run example.jst`. You should see output similar to the example below:

```shell
[🌵  jetstream] 2018-10-30 16:49:05: Version jetstream 1.0.2.dev0
[🌵  jetstream] 2018-10-30 16:49:05: Cmd args: /home/rrichholt/.local/bin/jetstream run example.jst 
[🌵  workflows] 2018-10-30 16:49:05: Building workflow...
[🌵  workflows] 2018-10-30 16:49:05: Workflow ready: <jetstream.Workflow Counter({'new': 2})>
[🌵     runner] 2018-10-30 16:49:05: Starting run: js01CV3P2XER2HTKGY5PP6JCJ37R
[🌵      local] 2018-10-30 16:49:05: LocalBackend initialized with 16 cpus
[🌵      tasks] 2018-10-30 16:49:05: <Task 51c2b5a7> has started
[🌵      local] 2018-10-30 16:49:05: Spawn: <Task 51c2b5a7>
[🌵      tasks] 2018-10-30 16:49:05: <Task e9dfdee4> has started
[🌵      local] 2018-10-30 16:49:05: Spawn: <Task e9dfdee4>
Tue Oct 30 16:49:05 MST 2018
[🌵      tasks] 2018-10-30 16:49:05: <Task 51c2b5a7> is complete
[🌵      local] 2018-10-30 16:49:05: Done: <Task 51c2b5a7>
dback-login1
[🌵      tasks] 2018-10-30 16:49:05: <Task e9dfdee4> is complete
[🌵      local] 2018-10-30 16:49:05: Done: <Task e9dfdee4>
[🌵     runner] 2018-10-30 16:49:06: Shutting down runner...
[🌵     runner] 2018-10-30 16:49:06: Run complete: js01CV3P2XER2HTKGY5PP6JCJ37R
```

## Dependency Directives in Templates

Dependencies between tasks can be specified with directives like `before` and `after`.
In this example, "task2" needs to run after "task1". Adding the `after` directive
will establish a dependency where "task2" _depends_on_ "task1". 

```yaml
# Example 1

  - name: task1
    cmd: hostname > info.txt

  - name: task2
    after: task1
    cmd: date >> info.txt

```

Or, equivalently:

```yaml
# Example 2

  - name: task1
    before: task2
    cmd: hostname > info.txt

  - name: task2
    cmd: date >> info.txt

```

Another way to set up dependencies is with the `input` and `output` directives. Example 
3 sets up the same workflow as Example 1, but using different directives.

```yaml
# Example 3

  - cmd: hostname > info.txt
    output: info.txt

  - cmd: date >> info.txt
    input: info.txt

```

> `output` directives can be declared without having a task with a matching `input` directive. But,
  an error will be raised if there are `input` directive with no tasks have a matching `output` directive.


All dependency directives (`before`, `after`, `input`, `output`) can be sequence of items. 
This example will setup a third task that waits for both setup tasks to complete 
before it executes:

```yaml
# Example 4

  - name: setup_task_1
    cmd: hostname > info.txt

  - name: setup_task_2
    after: setup_task_1
    cmd: date > info.txt

  - name: start
    after: 
      - setup_task_1
      - setup_task_2
    cmd: cat info.txt

```

> Note that there are multiple syntax options for a 
  [sequence in YAML](http://yaml.org/spec/1.2/spec.html#id2790320).


All dependency directives (`before`, `after`, `input`, `output`) can also be regex patterns.
Regex patterns can be a faster way of setting up dependencies in some cases. Here's the 
same workflow written with a pattern for the `after` directive instead of a sequence. 

```yaml
# Example 5 

  - name: setup_task_1
    cmd: hostname > info.txt

  - name: setup_task_2
    after: setup_task_1
    cmd: date >> info.txt

  - name: start
    after: setup_task_.*
    cmd: cat info.txt

```

> Sequences of patterns are fine too.


### Variables in templates

[Jinja2](http://jinja.pocoo.org) syntax can be used to add dynamic elements
to workflow templates. Dynamic elements are usually required to adapt tasks
to the configuration/input data for a pipeline. To keep it simple for now,
were going to build a toy pipeline: Sandwich Pipe. 

The goal of Sandwich Pipe is to build a sandwich. The commands are going to 
run a fake program called sandwich. And we're going to require a few 
configuration settings: `bread_type`, `meat_type`, `add_cheese`.

```yaml

  - name: get_bread
    cmd: sandwich get_bread --{{ bread_type }} > sandwich.txt

  - name: add_meat
    cmd: sandwich add_meat --{{ meat_type }} >> sandwich.txt
  
  {% if add_cheese %}
  # Sorry, we've only got one type of cheese
  - name: add_cheese
    cmd: sandwich add_cheese >> sandwich.txt
  {% endif %}

```

TODO: How to get configuration data into templates


This example only has two tasks written, but the loop will duplicate them, 
for a total of 6 tasks added to the workflow. Each "task2" will wait for its 
corresponding "task1" to complete. Also notice, "task3" has a regex pattern for 
its "after" directive. This will set up dependencies for every task matching the 
pattern. So, "task3" will wait for all "task2_A", "task2_B", and "task2_C" to
complete before starting.

```yaml

    {% for i in [A, B, C] %}

    - name: task1_{{ i }}
      cmd: hostname > info_{{ i }}.txt

    - name: task2_{{ i }}
      after: task1_{{ i }}
      cmd: date >> info_{{ i }}.txt

    {% endfor %}

    - name: task3
      after: task2_.*
      cmd: who
```

# Building Pipelines from Python Modules

TODO: This section is not added yet...

# Task Directives

Tasks are just a set of key-value properties called directives. If you're writing pipelines as a 
template, the task directives will look like: `directive: value`. If you're building tasks with 
the Python API, they would be key-word arguments to the function: `Task(directive=value)`. 
Some directives will be used differently depending on which backend you're running with: Local 
vs Slurm. The order of order of the directives inside a task is irrelevant.

Here is a list of the current task directives and how they're used: 

- `name` 

  Allows for linking dependencies to this task. Used
  to determine the filename for logs.

- `after` 

  This task will run after tasks matching each given value. 
  Supports sequences and Python regex patterns.
  
- `before`

  This task will run before tasks matching each given value. 
  Supports sequences and Python regex patterns.
  
- `input`:

  This task will run after tasks with output directives matching 
  each given value. Supports sequences and Python regex patterns.
    
- `output`: 

  This task will satisfy a matching input value requirement.
  Supports sequences and Python regex patterns.
    
- `cmd`: 

  Shell command executed with "/bin/bash"

- `exec`: 

  Python code that will execute in the runner process immediately 
  before sending the task to the backend for execution. This feature
  is not meant to be commonly used, and cmd should be the first
  choice.
    
- `stdin`: 

  cmd stdin will be connected to this value

- `stdout`: 

  cmd stdout will be connected to this value

- `stderr`: 

  cmd stderr will be connected to this value

- `cpus`: 

  LocalBackend - Will reserve local cpus when launching cmd
  SlurmBackend - Passed as "-c" when requesting the job allocation
      
- `mem`: 

  SlurmBackend - Passed as "--mem" when requesting job allocation

- `methods`: 

  Describe what this task does in plain language, later this can
  be used to generate a methods section for a project.
  
- `description`: 
  
  Describe this task for potential developers


# Installation

## Recommended: Install with Pip

> TGen users on Dback can load the latest version of Python with `module load python`.

This is a Python package and requires Python3. Installation guides for Mac/Windows/Linux are available from 
the [Hitchiker's Guide to Python][install_help] After Python3 is installed, you can install Jetstream with 
Pip (next step).


Install via HTTPS (requires entering Github username and password)

```shell
pip3 install --upgrade --user git+https://github.com/tgen/jetstream.git@master
```

Install via SSH (requires SSH keys to be configured in Github profile)

```shell
pip3 install --upgrade --user git+ssh://git@github.com/tgen/jetstream.git@master
```

##  Development install


Best practices for development are to clone the source repo, create a
virtual environment and install as an editable link with Pip. This will
cause (most) changes made in the source code to be testable without
resintalling. Here is an example command list:

```shell
   git clone https://github.com/tgen/jetstream.git
   virtualenv jetstream_venv
   source jetstream_venv/bin/activate
   pip3 install -e jetstream/
```

## Command not found after install

Pip installs package scripts and executables automatically. In some
cases you will need to update your [`PATH`][PATH] to include the directory where 
these scripts are stored.

- **MacOS**:

    Changes should be added to ``~/.bash_profile``.

    Installed system-wide: the bin directory will be ``/usr/local/bin/``

    Installed with ``--user``: the bin directory will be
    ``~/Library/Python/<pythonversion>/bin``

- **Linux**:

    Changes should be added to ``~/.bashrc``

    Installed system wide: the bin directory will be ``/usr/bin/``

    Installed with ``--user``: the bin directory will be ``~/.local/bin/``

[`See this post`][path_help] for more information.


# Tips

## Syntax Highlighting

Sublime text works well for writing templates. There is a syntax highlighter package for Ansible that 
will make it easier to read. First install the Ansible highlighter:

```
⌘-⇧-P
Package Control: Install Package
(Search for Ansible)
```

Then, configure the editor to use the "Ansible" syntax highlighter for `.jst` extensions:

`View -> Syntax -> Open all with current extension as ...`

Ansible is a system management utility that uses YAML/Jinja documents very similar to workflow 
templates. Documentation from Ansible covering the YAML template syntax can be found at: 
https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html


## Understanding YAML

YAML is a spec that allows common data structures in to be stored as unambiguous 
plain-text. YAML is the foundation of the workflow template syntax, so it's 
important to have a good grasp on the fundamentals. The full specification can 
be found [here](http://yaml.org), this section is just meant as a crash course 
introduction.

## Scalars

Scalars are single values of a specific type: Strings, Integers, Floats,
or Boolean. These basic types are present in almost every programming language
in some form or another. Yaml syntax provides a way to unambigously express the
type of a value in plain text: ``"true" vs true``  or ``"1234" vs 1234``

### Sequences

Are a way to represent an multiple items in ordered sets: list (Python), vector
(R), array (Javascript)

Python:

```python
chars = ['Philip J. Fry', 'Bender Bending Rodriguez', 'Leela Turanga']
```

Yaml:

```yaml
- Philip J. Fry
- Bender Bending Rodriguez
- Leela Turanga
```

### Mappings

Are a way to represent a set of paired key-values: dictionary (Python), named
vector (R), object (Javascript)

Python:

```python
pjf = {
    'first_name': 'Philip',
    'middle_name': 'J.',
    'last_name': 'Fry'
    'score': 42
}
```

Yaml:

```yaml
first_name: Philip
middle_name: J.
last_name: Fry
n_eyes: 2
```

### Or, combine of all of the above for maximum power

The values inside a sequence dont have to be scalars, they can be 
other sequences or mappings. Here is an example of a Sequence of 
Mappings written in Python. Notice here we're using lists and 
dictionaries.

```python
chars = [
    {
        'first_name': 'Philip',
        'middle_name': 'J.',
        'last_name': 'Fry',
        'n_eyes': 2
    },
    {
        'first_name': 'Bender',
        'middle_name': 'Bending',
        'last_name': 'Rodriguez',
        'n_eyes': 2
    },
    {
        'first_name': 'Leela',
        'middle_name': None,
        'last_name': 'Turanga',
        'n_eyes': 1
    }
]
```

The same structure can be created in plain text with YAML:

```yaml
- first_name: Philip
  middle_name: J.
  last_name: Fry
  n_eyes: 2

- first_name: Bender
  middle_name: Bending
  last_name: Rodriguez
  n_eyes: 2

- first_name: Turanga
  middle_name: null
  last_name: Leela
  n_eyes: 1
```

# Advanced Workflow Design

This page describes some of the more advanced features in Jinja2 that
can be used to write workflow templates.

## Filter Magic

Filters are used for applying transformations to the context data. Here
are some examples of filters that I've found really useful:

selectattr

Used inside of a `{% set ... %}` statement, this filter can be used to
create subsets of data. This example filters the read\_groups list for
data that is only Genome or
Exome:

``` sourceCode yaml
{% set dna_rgs = rgs | selectattr("General Library Type", "in", ["Genome", "Exome]) %}
```

## Macros

This is an example of how to develop workflow components that behave
more like the functions you would find in a traditional programming
language.

First lets make a basic macro. It works just like a function, allowing
for input arguments, and "returning" text to the template:

The first line starts the declaration, "gatk\_whatver" is the name of
the macro. "batch" is some iterable containing a set of intervals where
GATK should operate. "sample\_name" is the name of the sample were
running. When called, this macro will produce a GATK task that runs on
all the intervals in the batch.

``` sourceCode yaml
{% macro gatk_whatever(batch, sample_name) %}
- name: test_{{ sample_name }}
  cmd: bash
  stdin: |
    gatk Whatever \
    {% for i in batch %}
      -L {{ i }} \
    {% endfor %}
      --out {{ sample_name }}.out
{% endmacro %}
```

We can call the macro like this:

``` sourceCode yaml
{{ gatk_whatever(batch=[1,2,3], sample_name='foo') }}
```

## Handling whitespace with macros

Macros are a powerful way to create tasks, but handling whitespace can be a pain.
YAML is whitespace-sensitive and will throw errors if your task content does not
fit the proper indentation levels after rendering. This example shows a few ways
of calling macros that can be used to add or remove whitespace depending on your 
needs:

Here is a template that includes a simple macro and a few calls:
```
{% macro foo() %}
bar
  bar
    bar
  bar
bar
{% endmacro %}

Plain call:
{{ foo() }}

Indented call:
  {{ foo() }}

Call piped to indent(2):
{{ foo()|indent(2) }}

Indented call piped to indent(2):
  {{ foo()|indent(2) }}

Call with the first newline removed:
{{- foo() }}

Call with the last newline removed:
{{ foo() -}}

Last line of template
```

And here are the results of rendering: `jetstream render example.jst`
```

Plain call:
bar
  bar
    bar
  bar
bar


Indented call:
  bar
  bar
    bar
  bar
bar


Call piped to indent(2):
bar
    bar
      bar
    bar
  bar


Indented call piped to indent(2):
  bar
    bar
      bar
    bar
  bar


Call with the first newline removed:bar
  bar
    bar
  bar
bar


Call with the last newline removed:
bar
  bar
    bar
  bar
bar
Last line of template
```

## Macro-ception

Next is an example of a nested macro: a macro that calls another macro.
This is similar to a decorator in Python. This macro can be used to
group one list into smaller lists and then call some other function for
each batch. Arguments can be passed to the nested function, fn, by
including positional args in a list, or key-word args in a
dict.

``` sourceCode yaml
{% macro for_batch_n(intervals, fn, size=3, args=list(), kwargs=dict() ) %}

{% set batches = intervals | batch(size) | list %}
{% for batch in batches %}
{{ fn(batch, *args, **kwargs) }}
{% endfor %}

{%- endmacro %}

And here we call the batcher function, and tell it to call
"gatk_whatever" for each batch it produces.
```

``` sourceCode yaml
{{ for_batch_n(intervals, gatk_whatever, size=4, args=[sample_name,]) }}
```


[PATH]: http://www.linfo.org/path_env_var.html
[install_help]: http://docs.python-guide.org/en/latest/starting/installation/
[path_help]: https://stackoverflow.com/questions/35898734/pip-installs-packages-successfully-but-executables-not-found-from-command-line

