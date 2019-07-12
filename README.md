# Jetstream

Jetstream is a pipeline development framework written as a pure Python package 
and command line utility. It supports complex workflows modeled as directed-
acyclic graphs (DAGs), and execution on batch schedulers (Slurm).


Jetstream supports a variety of use cases to fit your individual needs:

- Run simple workflows like scripts for jobs that require just a few options

- Create projects that keep track of logs, configuration, and progress

- Or, design pipelines to solve complex problems with: customizable
  configuration files, check-pointing, and version control.


Simple workflows can be written as YAML documents that contain a set of _tasks_ 
to run. Tasks contain the commands to execute and attributes to control when 
and how they will be executed. For resuability, variables can be added to 
the document with the Jinja2

```yaml
- output: foo.txt
  cmd: echo hello world > foo.txt

- input: foo.txt
  cmd: cat foo.txt | say
 
```

And as the need for complexity increases, the expressive language of the 
templating system can be used to adapt the document to your input data:   

```yaml
{% for sample in samples %}
- name: haplotypecaller_{{ sample.name }}
  cmd: |
    activate gatk4

    gatk \
      -T HaplotypeCaller \
      -R "{{ reference_fasta }}" \
      -I "{{ sample.bam_path }}" \
      -o "{{ sample.name }}.raw.indels.snps.vcf"

{% endfor %}
```

# Usage

There are two ways to work with Jetstream: via the command-line utility or 
as a Python package. The majority of this guide will cover the command-line
interface which is the most common use case. See [installation](#installation) 
for help getting started.

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

- [**Templates**](#templates): Are text documents that can be processed and 
  run with the `jetstream` command-line utility.
  
- [**Pipelines**](#pipelines): Are templates that have been organized and 
  tagged with version information. They can be used with `jetstream pipelines` 


## Python package

```python
import jetstream
wf = jetstream.Workflow()
wf.new_task(name='task1', cmd='echo hello world')
wf.new_task(name='task2', after='task1', cmd='echo all done')
graph = wf.graph()
```

- [**Python Modules**](#building-pipelines-from-python-modules): Using the 
  functions and classes in the `jetstream` package to manually construct a
  workflow object with Python code.
  

# Templates

Jetstream models pipelines with _Tasks_ and _Workflows_:

- **Tasks** are the fundamental unit of a workflow. They describe _what to do_: 
  usually a shell command to execute. And _when to do it_: with flow directives 
  like `before`, `after`, or `input` that establish dependencies between tasks.

- **Workflows** are set of tasks along with information about their current 
  state

  Task "flow" directives (like input, output, before, and after) can be used to 
  identify dependencies between tasks, and construct a DAG model workflow for a 
  project. Most users will not need to interact directly with workflow files, 
  they will be handled by the command-line utility.

Templates are a way to create a set of reusable _tasks_ that can run with 
multiple sets of input data. You can think of templates like a declarative 
scripting language for pipelines:

```
# example.jst

- name: say_hello
  cmd: say "hello {{ name }}"
  
```

They can be run with the command-line tool:

```
$ jetstream run example.jst -c name bender
```

**Templates** are text files that describe the tasks in a pipeline. They 
(usually) have variable components that will be filled in with project-specific 
data at runtime. Once the template is combined with project-specific 
configuration data, the tasks will be a "rendered" and saved as a workflow.

Each task starts with a hyphen, `-`, and then task directives are given 
as `directive: value` lines. This is YAML syntax for describing a sequence of 
mappings. If you're unfamiliar with YAML, see the 
[introduction in docs](docs/yaml_help.md).

Here is an example template that describes two tasks to complete:

```yaml

- name: task1
  cmd: hostname

- name: task2
  cmd: date

```

In this case the tasks are relatively simple: they have a `name` that allows us 
to link dependencies, and a `cmd` that will be executed. There are no variables 
or flow directives yet, we'll add those in the next sections. Since there are 
no dependencies declared, the runner may execute both tasks in parallel if the 
computing resources are available


## Dependencies between Tasks in templates

Dependencies between tasks can be specified with "flow" directives: `before`, 
`after`, `input` and `output`. In this example, "task2" needs to run after 
"task1". Adding the `after` directive will establish a dependency where "task2" 
_depends_on_ "task1". 

```yaml

- name: task1
  cmd: hostname > info.txt

- name: task2
  after: task1
  cmd: date >> info.txt

```

Another way to set up dependencies is with the `input` and `output` directives. 
This example sets up the same workflow as above, but using different directives.

```yaml

- cmd: hostname > info.txt
  output: info.txt

- cmd: date >> info.txt
  input: info.txt

```

> `output` directives can be declared without having a task with a matching 
  `input` directive. But, an error will be raised if there are `input` 
  directives with no matching `output` directives.


All dependency directives (`before`, `after`, `input`, `output`, and `*-re` 
variants) can be sequences. This example will setup a third task that waits for 
both setup tasks to complete before it executes:

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


Some dependency directives (`before-re`, `after-re`, `input-re`) are treated as 
regex patterns. Regex patterns can be an easier way of setting up dependencies 
when there are many tasks you need to link. Here's the same workflow written 
with a pattern for the `after-re` directive instead of a sequence. 

Note that patterns and input/output directives will increase the time it takes
to calculate the workflow graph. For most workflows this is negligible but may
become a concern for very large workflows.

```yaml

- name: setup_task_1
  cmd: hostname > info.txt

- name: setup_task_2
  after: setup_task_1
  cmd: date >> info.txt

- name: start
  after-re: setup_task_.*
  cmd: cat info.txt

```

similarily input/outputs can be used:

```yaml

- name: setup_task_1
  output: info1.txt
  cmd: hostname > info1.txt

- name: setup_task_2
  after: setup_task_1
  output: info2.txt
  cmd: date >> info.txt

- name: start
  input-re: ".*\.txt"
  cmd: cat info.txt

```


## Variable data in templates

[Jinja2](http://jinja.pocoo.org) can be used to add variables to workflow 
templates. Prior to loading tasks and connecting dependencies, Jetstream will 
render templates with Jinja2. During the render, variables in the template are 
replaced with actual data given as command arguments, config files, or saved 
in the project and pipeline indicies. Templating is a pattern used widely in 
web development, but other examples of using Jinja and YAML together to create
dynamic structured documents can be found in tools like Ansible or SaltStack. I 
recommend exploring tutorials on Jinja2 in addition to reading this section.

- Add variables to templates with the double-curly-bracket syntax:

  ```
  - name: say_hello
    cmd: say "hello {{ name }}"

  ```

- Use logical expressions with curly-bracket-parenthesis syntax:

  ```
  # comments are allowed throughout these documents following either yaml or
  # jinja2 comment syntax

  {% for name in names %}
  - name: say_hello
    cmd: say "hello {{ name }}"

  {% endfor %}
  ```
  
- And much, much more... see 
  [designer documentation for details](http://jinja.pocoo.org/docs/2.10/templates/)

When the template is rendered, data used for filling in variables is pulled 
from several sources. Each is explained below in further detail:
 
1) _command-line arguments_: `-c/--config` or `-C/--config-file`
2) data stored in the _project_ index (if using a project): 
   `<project>/jetstream/project.yaml`
3) data stored in the _pipeline_ manifest (if using a pipeline): 
   `<pipeline>/pipeline.yaml`


After these sources are loaded, they're collapsed into a single config object
(a dictionary) that is used by Jinja2 as the _context_ for rendering the 
template.
 
#### Template Data from Command-line Arguments

In the example below, the variable `{{ name }}` is what we want to replace, so 
we need to pass in config data with key `name`. To pass a single variable with 
command line arguments, use the `-c/--config` options. It takes two arguments:
the first is the key, the second is the value.  The key can optionally include 
the type of the variable being passed with the syntax: `type:key value`. If
type is not given, the data will be loaded as a string. This syntax supports a 
wide range of datatypes, including entire files. 

These arguments are also allowed when creating projects with `jetstream init` 

You can test this process with the `jetstream render` command. The resulting
template will be printed to stdout:

```yaml
# example1.jst
- name: say_hello
  cmd: say "hello {{ name }}"

```

```bash
$ jetstream render example1.jst -c name bender
# example1.jst

- name: say_hello
  cmd: say "hello bender"
  
```

Batches of config data can also be loaded from files. There is a dynamic file
loader type that will handle json, yaml, and many tabular text file formats. It
will determine the file type based on the extension of the path. Here the
example template has been modified to accept a set of names. We can load that
set of names from a json file:

`example2.jst`
```yaml
# example2.jst
{% for name in names %}
- name: say_hello
  cmd: say "hello {{ name }}"

{% endfor %}
```

`config.json`
```
{"names": ["Philip J. Fry", "Bender Bending Rodriguez", "Leela Turanga"]}
```

To run, use `-C/--config-file` to load the entire file of variables.

```bash
$ jetstream render example2.jst -C config.json
# example2.jst
- name: say_hello
  cmd: say "hello Philip J. Fry"

- name: say_hello
  cmd: say "hello Bender Bending Rodriguez"

- name: say_hello
  cmd: say "hello Leela Turanga"
 
```

Both methods can be used together as well, but loading a `--config-file` 
overwrite any other config data loaded from previous arguments, so it should
usually be given first.


#### Template Data saved in Projects

_Projects_ are an optional but very helpful feature in Jetstream. A project is 
a directory that contains a jetstream folder and `project.yaml` (this folder 
will be referred to as the project _index_). You can create these directories 
with the `jetstream init` command.

> Config variables can be used during project init. Any data will be saved into
  the project index so that it is available when running templates or 
  pipelines on that project in the future.
  
When running many jetstream commands (`project`, `tasks`, `run`, `pipelines`, 
etc.) projects help jetstream to organize task data, logs, and store workflow
progress. Projects also serve to store configuration data that will be used
when running templates and pipelines. In the project index (`jetstream` folder)
there is a `project.yaml` file that contains info about the project when it 
was created. This file is always included as a config data source when 
rendering templates or running pipelines.


#### Template Data saved in Pipelines

_Pipelines_ may specify additional data that is available when rendering their 
templates. Pipelines should include a `pipeline.yaml` file. Inside this file
the required field `__pipeline__` contains information about the pipeline, but
any additional fields can be specified. 

An example use case for this would be when if a pipeline had a set of possible
options to choose from. They could be stored in the `pipeline.yaml`, and remove
the amound of config data needed to be included for each invocation of the 
pipeline. Here's an example:

`pipeline.yaml`
```yaml
__pipeline__:
  name: example_pipe
  ...
reference_file: gs://bucket/path/to/reference.file
run_modes:
  default:
    threads: 4
    index_uri: gs://bucket/path/to/basic.index
  faster:
    threads: 16
    index_uri: gs://bucket/path/to/faster.index
    
``` 


# Tasks

Tasks are just a set of key-value properties called directives. If you're 
writing pipelines as a template, the task directives will look like this: 
`directive: value`. Some directives will be used differently depending on 
which backend you're running with: Local vs Slurm. The order of order of 
the directives inside a task is irrelevant. 

Any additional directives can be added as well, they will not affect the runner
or execution of tasks. These will be stored in the workflow and could be useful 
for downstream meta-analysis of the projects themselves (eg. performance 
comparisons, or documenting the tasks)

Here is a list of the common task directives and how they're used: 

---

"Core" directives

- `name` 

  A unique identifier for the task. If this is absent, it will be assigned 
  based on a hash of the task content. Allows for linking dependencies to this 
  task. Used to determine the filename for logs.
  
- `cmd`: 

  Command to be executed by the backend, should be valid 
  [Bash](https://www.gnu.org/software/bash/).
  This is where the main "work" of the task should exist.

---

"Flow" directives 

- `after` 

  This task will run after tasks named with each value. Supports sequences.
  
- `after-re`
  
  This task will run after tasks matching each given regex pattern. 
  Supports sequences.
  
- `before`

  This task will run before tasks named with each value. Supports sequences.

- `before-re`
  
  This task will run before tasks matching each given regex pattern. 
  Supports sequences.
  
- `input`:

  This task will run after tasks with output directives matching 
  each given value. Supports sequences.

- `input-re`:

  This task will run after tasks with output directives matching 
  each given regex pattern. Supports sequences.
  
- `output`: 

  This task will satisfy a matching input value requirement.
  Supports sequences.

---

"Execution" directives

- `exec`: 

  Python code that will execute in the runner process immediately 
  before sending the task to the backend. This feature can be used
  to modify the workflow (add tasks) while it's being run. Two local
  variables are added during execution: `task` and `runner`. The runner
  is important because it contains the current workflow. Any
  errors during execution will halt the runner immediately.
  
  Note: the workflow graph will always be recalculated after any exec 
  directive runs, and most work can be done with `cmd` directives. 
      
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


# Pipelines

Pipelines are Jetstream templates that have been documented with version
information and added to a jetstream pipelines directory. The 
`jetstream pipelines` command allows pipelines to be referenced by name and 
automatically includes any pipeline scripts and variables in the run.

To create a pipeline: 
    
- Add the template file(s) to a directory that is in your pipelines 
  searchpath. The default searchpath is your user home directoy, but it 
  can be changed in the application settings (see `jetstream settings -h`) 

- Create a `pipeline.yaml` file in the directory and be sure to include the 
  required fields. Here is an example: 
  
`pipeline.yaml` 
```yaml
# Required info about the pipeline
__pipeline__:
  name: # Alphanumeric pipeline name
  version: # version number should be parsable with Python looseversion
  main: # the path (relative to this file) to the template file
  bin:  # *Optional* this directory will be prepended to env variable PATH 

# *Optional* any other fields in this file
foo: 42
bar: 24
``` 

# Installation

## Recommended: Install with Pip

> TGen users on Dback can load the latest version of Python with 
  `module load python`.

This is a Python package and requires Python3. Installation guides for 
Mac/Windows/Linux are available from the 
[Hitchiker's Guide to Python][install_help] After Python3 is installed, you can 
install Jetstream with Pip (next step).


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
cases you will need to update your [`PATH`][PATH] to include the directory 
where these scripts are stored.

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

[PATH]: http://www.linfo.org/path_env_var.html
[install_help]: http://docs.python-guide.org/en/latest/starting/installation/
[path_help]: https://stackoverflow.com/questions/35898734/pip-installs-packages-successfully-but-executables-not-found-from-command-line

