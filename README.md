# Jetstream

Jetstream is a pipeline development framework written as a pure Python package 
and command line utility. It supports complex workflows modeled as directed-
acyclic graphs (DAGs), and execution on batch schedulers (Slurm).

```
- cmd: echo Hello, World
```

There are two ways to work with Jetstream: via the command-line utility or 
as a Python package. The majority of this guide will cover the command-line
interface which is the most common use case. See [installation](#installation) 
for help getting started.

---

# Use Jetstream

Jetstream supports a variety of use cases to fit your individual needs:

- Run simple workflows like scripts for jobs that require just a few options

- Create projects that keep track of logs, configuration, and progress

- Or, design pipelines to solve complex problems with: customizable
  configuration files, check-pointing, and version control.


## Simple workflows

Workflows are organized with YAML and contain a set of commands to run along
properties that control when and how they will be executed. Here is an example


```yaml
- name: first_task
  cmd: echo Hello World

  
- name: second_task
  cmd: echo 

```

## Command-line

After installing with Pip, the command line utility can be launched with
`jetstream`, and help can be accessed with the `-h/--help` options. If 
the command is not found, see the detailed [installation](#installation) 
help below.

```shell
$ jetstream -h
```

The command line utility includes tools for managing projects, building 
pipelines, and executing runs. The most common use case is designing a 
reusable template that can process data for many different projects:

- [**Templates**](#building-pipelines-from-templates): Are text documents
  that can be processed and run with the `jetstream` command-line utility.
  

## Python package

```python
>>> import jetstream
>>> wf = jetstream.Workflow()
>>> wf.new_task(name='task1', cmd='echo hello world')
>>> wf.new_task(name='task2', after='task1', cmd='echo all done')
>>> graph = wf.graph()
```

- [**Python Modules**](#building-pipelines-from-python-modules): Using the 
  functions and classes in the `jetstream` package to manually construct a
  workflow object with Python code.
  

# Templates

Jetstream models pipelines with _Tasks_ and _Workflows_:

- **Tasks** are the fundamental unit of a workflow. They describe _what to do_: 
  usually a shell command to execute. And _when to do it_: with flow directives 
  like `before`, `after`, or `input` that establish dependencies between tasks.

- **Workflows** are set of tasks along with information about their current state

  Task "flow" directives (like input, output, before, and after) can be used to identify 
  dependencies between tasks, and construct a DAG model workflow for a project. Most users
  will not need to interact directly with workflows, they will be handled by the command-line
  utility.

Templates are a way to create a set of reusable _tasks_ and run them on multiple projects.
Think of them like a declarative scripting language for pipelines:

```
# example.jst

- name: say_hello
  cmd: say "hello {{ name }}"
  
```

They can be run with the command-line tool:

```
$ jetstream run example.jst -- --str:name bender
```

**Templates** are text files that describe the tasks in a pipeline. They (usually) have 
variable components that will be filled in with project-specific data at runtime. 
Once the template is combined with project-specific configuration data, the tasks will be a 
"rendered" and ready to run.

Each task starts with a hyphen, `-`, and then task directives are given as `directive: value` 
lines. This is YAML syntax for describing a sequence of mappings. If you're unfamiliar with YAML, 
see the [introduction below](#understanding-yaml).

Here is an example template that describes two tasks to complete:

```yaml

- name: task1
  cmd: hostname

- name: task2
  cmd: date

```

In this case the tasks are relatively simple: they have a `name` that allows us to link 
dependencies, and a `cmd` that will be executed. There are no variables or flow directives
yet, we'll add those in the next sections. 


## Dependencies between Tasks in templates

Dependencies between tasks can be specified with "flow" directives: `before`, `after`, 
`input` and `output`. In this example, "task2" needs to run after "task1". Adding the
`after` directive will establish a dependency where "task2" _depends_on_ "task1". 

```yaml

- name: task1
  cmd: hostname > info.txt

- name: task2
  after: task1
  cmd: date >> info.txt

```

Another way to set up dependencies is with the `input` and `output` directives. This example
sets up the same workflow as above, but using different directives.

```yaml

- cmd: hostname > info.txt
  output: info.txt

- cmd: date >> info.txt
  input: info.txt

```

> `output` directives can be declared without having a task with a matching `input` directive. But,
  an error will be raised if there are `input` directives with no matching `output` directives.


All dependency directives (`before`, `after`, `input`, `output`) can be sequences. 
This example will setup a third task that waits for both setup tasks to complete 
before it executes:

```yaml

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

- name: setup_task_1
  cmd: hostname > info.txt

- name: setup_task_2
  after: setup_task_1
  cmd: date >> info.txt

- name: start
  after: 
    - re: setup_task_.*
  cmd: cat info.txt

```

> Sequences of patterns are fine too.


## Variables in templates

[Jinja2](http://jinja.pocoo.org) can be used to add variables to workflow templates.
Prior to loading tasks and connecting dependencies, Jetstream will render a template
with Jinja2. During the render, variables in the template are replaced with actual 
data given as command arguments, files, or project config files. This pattern of is widely
used in web development, but tools like Ansible and SaltStack also use Jinja+YAML to 
create dynamic, structured documents. I recommend exploring tutorials on Jinja2 in addition 
to reading this section.

- Add variables to templates with the double-curly-bracket syntax:

  ```
  # example1.jst

  - name: say_hello
    cmd: say "hello {{ name }}"

  ```

- Use expressions with curly-bracket-parenthesis syntax:

  ```
  # example2.jst

  {% for name in names %}
  - name: say_hello
    cmd: say "hello {{ name }}"

  {% endfor %}
  ```

When the template is run, variables can be filled in with _command-line arguments_, 
configuration data stored in the _project_, or user/pipeline _constants_.

#### Template Data from Command-line Arguments (Kvargs)

Any arguments following `--` will be treated as template variable data args. 
They should follow the syntax: `--<type>:<key> <value>` In the example above,
the variable `{{ name }}` is what we want to fill in, so the key `name` was
used. The type is `str`, handled by Python `str()`. And the value `bender` was
passed to that function. 

This syntax supports a wide range of datatypes, including entire files. More 
details about the loaders can be found in the kvargs module. 

You can test this process with the `jetstream render` command. The resulting
template will be printed to stdout:

```
$ jetstream render example1.jst -- --str:name "bender"
# example1.jst

- name: say_hello
  cmd: say "hello bender"
  
```

These arguments are also allowed when creating projects with `jetstream init`

#### Template Data from Project Config Files

_Projects_ are an important optional feature in Jetstream. A project is a directory
that contains a jetstream folder and `project.yaml`. You can create these directories
with the `jetstream init` command. 

Any arguments following `--` are treated as kvargs and will be stored in the project
config file. Running/rendering/building any templates in the project will load 
variables from the config file.

# Building Pipelines from Python Modules

TODO: This section is not added yet...


# Task Directives

Tasks are just a set of key-value properties called directives. If you're writing pipelines as a 
template, the task directives will look like this: `directive: value`. Some directives will be 
used differently depending on which backend you're running with: Local vs Slurm. The order of 
order of the directives inside a task is irrelevant. Here is a list of the current task directives 
and how they're used: 

---

"Core" directives

- `name` 

  A unique identifier for the task. If this is absent, it will be assigned based on a hash of the 
  task content. Allows for linking dependencies to this task. Used to determine the filename 
  for logs.
  
- `cmd`: 

  Command to be executed by the backend, should be valid 
  [Bash](https://www.gnu.org/software/bash/).
  This is where the main "work" of the task should exist.

---

"Flow" directives 

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

---

"Execution" directives

- `exec`: 

  Python code that will execute in the runner process immediately 
  before sending the task to the backend. This feature can be used
  to modify the workflow (add tasks) while it's being run.
    
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

---

"Description" directives

- `tags`: 
  
  Sequence of short descriptive tags
  
- `methods`: 

  Describe what this task does in plain language, later this can
  be used to generate a methods section for a project.
  
- `description`: 
  
  Describe this task for potential users


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

---

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

## Sequences

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

## Mappings

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

## Or, combine of all of the above for maximum power

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

Filters are used for applying transformations to the template data.
There are loads of filters available for use in templates, find info
here: http://jinja.pocoo.org/docs/2.10/templates/#builtin-filters

Here are some examples of filters that I've found really useful:

- `join`

  Join together a sequence of items:

  ```yaml
  - cmd: cat {{ filepaths|join(' ') }} > merged.txt
  ```

- `lower`/`upper`
  
  Make your variables case-insensitive with these handy guys:
  
  ```yaml
  {% if textfile.compression_type|lower == 'gz' %}
  gunzip -c {{ textfile.path }} > temp.txt
  {% endif %}
  ```
  
- `selectattr`

  Used inside of a `{% set ... %}` statement, this filter can be used to
  create subsets of data. This example filters the rgs list for data that 
  is only Genome or Exome:

  ```yaml
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

```yaml
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

