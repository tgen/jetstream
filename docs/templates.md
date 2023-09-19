# Templates

The building blocks of a Jetstream pipeline are [templates](docs/templates.md). 
Templates are scripts that describe the steps required for a pipeline. Some
pipelines can be a single template file. Other pipelines will require several 
templates and supporting data files (see [pipelines](docs/pipelines.md)).

# Table of contents

- [Templates](#templates)
- [Table of contents](#table-of-contents)
- [Syntax](#syntax)
- [Dependencies between _Tasks_ in templates](#dependencies-between--tasks--in-templates)
- [Variables and Logic](#variables-and-logic)
  * [Additional globals and filters](#additional-globals-and-filters)
    + [Globals](#globals)
    + [Filters](#filters)
- [Template rendering data](#template-rendering-data)
  * [Template Data from Command-line Arguments](#template-data-from-command-line-arguments)
  * [Template Data saved in Projects](#template-data-saved-in-projects)
  * [Template Data saved in Pipelines](#template-data-saved-in-pipelines)
- [Modularity](#modularity)

# Syntax 

Templates describe a set of [tasks](tasks.md) that can be run for multiple projects. 
You can think of templates like a declarative scripting language for building 
pipelines. They usually have variable components that will be filled-in with 
project data at runtime. Once the variables are filled-in, the tasks will be 
*rendered*, or finalized commands ready to execute.

Here is an example workflow template with one task to complete:

```yaml
- name: hello_world
  cmd: echo "{{ greeting }} world"
  
```

Here is a template with two tasks:

```yaml
- name: task1
  cmd: echo "{{ greeting }} world"

- name: task2
  after: task1
  cmd: echo "It is $(date)"

```

And here is a template with 100 tasks:

```yaml
{% for i in range(99) %}
- name: task1
  cmd: echo "{{ greeting }} world number {{ i }}"

{% endfor %}

- name: task2
  after: task1
  cmd: echo "It is $(date)"

```

The foundation of template syntax is [YAML](docs/yaml_help.md). Each task 
starts with a hyphen, `-`, and then a set of task directives are 
given as `directive: value` lines. This is YAML syntax for describing a 
sequence of mappings.

Directives are the instructions that the runner will use to determine _when_
and _how_ to execute the task. See [Tasks](tasks.md) for more information.

Variables and logic can be included to make templates dynamic. In the example
above, the `{{ greeting }}` will be filled in with data provided by the user
when the template is run with `jetstream run`. 
[Variables and logic](#Variable-and-logic) provide endless potential 
for modeling complex workflows.

Templates can be run with the command-line tool:

```
$ jetstream run example.jst -c name bender
```


# Dependencies between _Tasks_ in templates

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


# Variables and Logic

> Details for variables and logic syntax can be found here 
  [designer documentation for details](http://jinja.pocoo.org/docs/latest/templates/)

[Jinja2](http://jinja.pocoo.org) can be used to add variables and logic to
workflow templates. Prior to loading tasks and connecting dependencies,
Jetstream will render templates with Jinja2. During the render, variables in
the template are replaced with actual data given as command arguments, config
files, or saved in the project and pipelines. Templating is a pattern used 
widely in web development, but other examples of using Jinja and YAML
together to create dynamic structured documents can be found in tools like
Ansible or SaltStack. 

Examples:

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
  [designer documentation for details](http://jinja.pocoo.org/docs/latest/templates/)

## Additional globals and filters

In addition to the included global functions and filters included with Jinja2, 
several other tools have been added with Jetstream and can be used inside templates:

### Globals

- `raise`: Raise an error while rendering the template
  - Example: `{{ if foo < 42 }}{{ raise('foo should be at least 42') }}{{ endif }}`

- `log`: Log messages to the Jetstream logger while template renders
  - Example: `{{ log('Foo is {}'.format(foo), level='CRITICAL') }}`

- `env`: Returns environment variable value
  - Example: `echo foo is {{ getenv('FOO') }}`

- `getenv`: Returns environment variable value, this will return None if value 
  is not set whereas `env` will raise an error. A different fallback value can
  be given as the second argument.
  - Example: `echo foo is {{ getenv('FOO', None) }}`

- `setenv`: Sets an environment variable when the template is rendered
  - Example: `{{ setenv('FOO', '42') }}`


### Filters

- `fromjson`: Parse a json string as an object

- `basename`: Returns the basename of a path

- `dirname`: Returns the directory name of a path

- `urlparse`: Parse a url string as an object

- `sha256`: Returns sha256 hexdigest for a string

- `md5`: Returns md5sum of a file defined with a path
  - Example: `{{ required_scripts.some_script.path | md5 }}`

- `assignbin`: Returns the 0-based bin the value falls in. 
  - The default bin edges are 0 to infinity, meaning this will return 0 if the bin edges are not defined. 
  - Returns -1 if the input value is out of bounds.
  - Any value landing on an edge will floor to lower bin. 
  - This also accepts a list of labels such that: `{{ assignbin(5,[0,2,4,6],['low','med','high']) }}` returns 'high'. Moreover, `{{ assignbin(4,[0,2,4,6],['low','med','high']) }}` would return 'med'.

# Template rendering data

When the template is rendered, data is pulled from several sources. Each is
explained below in further detail (highest priority first):
 
1) _command-line arguments_: `-c/--config` or `-C/--config-file`

2) data stored in the _project_ index (if using a project): 
   `<project>/jetstream/project.yaml`

3) data stored in the _pipeline_ manifest (if using a pipeline): 
   `<pipeline>/pipeline.yaml`


After these sources are loaded, they're collapsed into a single config object
(a dictionary) that is used by Jinja2 as the _context_ for rendering the 
template. Higher-priority data sources will overwrite other sources.


## Template Data from Command-line Arguments

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

Batches of config data can also be loaded from files. Here the example template 
has been modified to accept a set of names. We can load that set of names from 
a json file:

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

There is a dynamic file loader type that will handle json, yaml, and many 
tabular text file formats. It will determine the file type based on the 
extension of the path. But, this can be overridden with the 
`--config-file-type` option. See `jetstream render -h` help for a list of 
supported file types for your configuration.


## Template Data saved in Projects

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


## Template Data saved in Pipelines

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

# Modularity

Templates can be modularized, or divided into smaller pieces, to improve 
organization and reusablility. There are a few ways to modularize code: 
`include`, `extends`, and `macros` are some options. Full details can be found 
in the Jinja2 documentation on 
[template inheritance.](https://jinja.palletsprojects.com/en/2.10.x/templates/#template-inheritance)
Here is an example using the `{% include %}` statement to use code from another
template file:

```yaml
# This template includes code from the next template

{% for sample in samples %}

- name: download_{{ sample }}
  cmd: wget $DOWNLOAD_URL_ROOT/{{ sample }}.gz

{% include 'process_sample.jst' with context %}

{% endfor %}

- name: finalize
  after: .*
  cmd: echo All done

```


```yaml
# This template code is used in the template above

- name: decompress_{{ sample }}
  cmd: gunzip {{ sample }}.gz


- name: transform_{{ sample }}
  after: decompress_{{ sample }}
  cmd: sed -i 's/monday/friday/' {{ sample }}.txt


- name: compress_{{ sample }}
  after: transform_{{ sample }}
  cmd: gzip {{ sample }}.txt

```
