# Pipelines

Pipelines are a method for packaging Jetstream templates in order to stay  
organized and improve usability. Pipelines encourage programming best-practices 
like documenation, modularization, and versioning.

Jetstream pipelines are directories with a pipeline manifest file. To be 
considered a pipeline, the directory needs to include the workflow templates and
also a [manifest file](#pipeline-manifest-files): `pipeline.yaml`. The directory
can also contain any supplementary files like documentation, reference data,
binaries, scripts, etc. Pipelines located in the
[searchpath](#pipeline-searchpath-lookup), can be run by name with the
`jetstream pipelines <name>` command. Pipelines can also be used with many other
Jestream commands via the `-p/--pipeline` option. Multiple versions of the same
pipeline can also be installed and referenced with: `jetstream pipelines
<name>@<version>`.

Here is an example of the directory tree for a `hello_world` pipeline. The 
contents of each file are also shown below:

```
~/hello_world/
  ├─ hello_world.jst
  ├─ pipeline.yaml
  └─ ... any additional files can be included here

```

`hello_world.jst`
```

- name: hello_world
  cmd: echo "{{ greeting }} world"
  
```

`pipeline.yaml`

```yaml
__pipeline__:
  name: hello_world
  version: 0.1
  main: hello_world.jst
```


## Pipeline searchpath lookup

When using the `jetstream pipelines` command, a search is performed for any
pipeline matching the name and optionally version number. If version number is
omitted, the latest version will be used. Jetstream looks for pipelines in 
any directory included in the pipelines *searchpath*. The default location is 
the user home directory, but this value can be changed in the user config file 
(see `jestream settings`). Pipelines can also be nested inside any other 
pipeline directory. The *searchpath* behaves much like the \*nix `PATH`
variable: multiple locations can be specified by separating with colons 
`<first path>:<second path>`, the first pipeline matching the search criteria
will be returned.


## Pipeline Manifest files

Pipeline manifest files must be saved in the top directory of the pipeline and
be named `pipeline.yaml`. Manifests must be a YAML formatted document with at
least a `__pipeline__` section.

The `__pipeline__` section is used by Jetstream to parse metadata about the 
pipeline. When running pipelines with the `jestream pipelines` command, version
information can be inspected to choose the latest (or specific) version to run.

The following values are required in the `__pipeline__` section:

- name: alpha-numeric and underscores allowed
- version: version parsing is loose, but a version must be present
- main: the name of the main template file that should be rendered when running
  the pipeline

Any additional items present in the pipeline manifest will be available when 
rendering templates. These values can be overridden by values added to a 
specific project, or by values provided via command-line arguments. These extra
fields are useful for providing defaults for pipeline behavior that can
be overridden via command-line arguments. Here is an example where a pipeline
may include default values for a variable by using the pipeline manifest:



`hello_world.jst`
```
- name: hello_world
  cmd: echo "{{ greeting }} world"
  
```

`pipeline.yaml`

```yaml
__pipeline__:
  name: hello_world
  version: 0.1
  main: hello_world.jst
greeting: hello
```

When running the pipeline, this value could be overridden:

```
$ jetstream pipelines hello_world -c greeting hola
```


## Pipeline environment variables

When running pipelines, several environment variables will be exported:

- `JS_PIPELINE_PATH`: The path to the pipeline
- `JS_PIPELINE_NAME`: The name of the pipeline
- `JS_PIPELINE_VERSION`: The version of the pipeline

These can be useful for introspection, logging, or using extra files packaged
with the pipeline. For example, a task command could reference data files
included in the pipeline:

```
- name: cat_pipeline_data_file
  cmd: |
    echo "Running ${JS_PIPELINE_NAME} version: ${JS_PIPELINE_VERSION}"
    cat "${JS_PIPELINE_PATH}/data/file.txt"

```

It's also a great way to keep pipeline code short by moving long commands into 
external scripts:

```
- name: run_pipeline_script
  cmd: |
    python "${JS_PIPELINE_PATH}/scripts/foo-v1.py"

```

If a pipeline includes a `bin` directory, it will be prepended to the `PATH` 
environment variable before running the pipeline. This allows any executables
to be referenced by name
