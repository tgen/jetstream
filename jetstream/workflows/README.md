# Workflows

The workflows module contains the `Workflow()` class and supporting functions.

# Workflow()

The `Workflow()` class models workflows as a directed-acyclic graph where nodes
are tasks to complete, and edges represent dependencies between those tasks. It
has methods for building workflows (add_node, add_dependency) in addition to
the methods required for executing a workflow (__next__, fail,
complete, etc.). Workflows can also be loaded from node-link data.


Workflows can also be described in a template syntax. Workflow templates can
contain variables to be filled in with project-specific data prior to loading.

Template + data --render--> Workflow()

# Details

Templates are an array of node described in YAML format:

```yaml

- id: Index Bams
  cmd: ['samtools', 'index', 'sampleA.bam']

```

On top of this, Jinja templating is used to add dynamic elements.

```yaml

{% for sample in project.config.samples %}
- id: Index Bam {{ sample.name }}
  cmd: ['samtools', 'index', '{{ sample.name}}.bam'
{% endfor %}

```


# Upfront workflow rendering

Workflows are an immutable data structure. Rendering a template is a dynamic
procedure that responds to input data and template directives to generate a
workflow. But, the resulting workflow is a final, and complete, description of
the tasks required, and the order of execution.


What about feedback? Conditionals?

The workflow is an immutable network graph defined prior to runtime.
It cannot be modified by events that occur during runtime. If a node exists
in a workflow, the runner will always launch it. The only exceptions to this
rule occur if 1) the runner exits before reaching that node, or 2) one of the
node's dependencies has failed.

Cases where feedback may be necessary:

Dynamic chunking of data

Some input data needs to be split into n chunks where n is determined
during workflow runtime. Each chunk then needs to be treated as an
individual task in the workflow by downstream tasks.

Note that this is not a problem if n can be determined prior to runtime,
or if your application can handle the chunking internally.


# Parsing/Serializing

json, yaml, etc..
