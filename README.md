# Jetstream

Jetstream is a pipeline development framework written as a pure Python package 
and command line utility. It supports complex workflows modeled as directed-
acyclic graphs (DAGs), and execution on batch schedulers (Slurm).


Jetstream supports a variety of use cases to fit your individual needs:

- Run simple workflows like scripts for jobs that require just a few options

- Create projects that keep track of logs, configuration, and progress

- Or, design pipelines to solve complex problems with: customizable
  configuration files, check-pointing, and version control.

# Table of contents

# Usage

Jetstream can be installed with pip. See [installation](#installation) for more
info:

```shell
pip install git+https://github.com/tgen/jetstream.git@master
```

After installing, there are two ways to work with Jetstream: 1) via the 
command-line utility or 2) as a Python package. The majority of this guide 
will cover the command-line interface which is the most common use case. 

After installing with Pip, the command line utility can be launched with
`jetstream`, and help can be accessed with the `-h/--help` options. If 
the command is not found, see the detailed [installation](#installation) 
help below.

```shell
$ jetstream -h
```


# Introduction

Jetstream pipelines are built by creating templates, then running those templates
by providing input data and a destination project. Simple pipelines may just need
to be a single template file, more complex use cases may have several template files
stored together in a directory with a `pipelines.yaml` manifest file. Detailed 
documentation for templates and pipelines can be found in the docs.

Templates are a way to create a set of reusable _tasks_ that can run with multiple 
sets of input data. You can think of templates like a declarative scripting language 
for building pipelines. Here is an example template file:

```
# save as example.jst

- name: hello_world
  cmd: echo "hello {{ name }}"
  
```

It can be run with the command-line tool:

```
$ jetstream run example.jst -c name world
```


# Intro Vignette

As a basic introduction, we'll create and run a pipeline that performs somatic 
variant calling on genome sequencing data. Jetstream template files are the 
building blocks for pipelines. They are simple text documents that containin a 
set of _tasks_ to run. It may help to think of templates like scripts that 
Jetstream can interpret and execute. 

* To get started, copy the following code to a new file. We'll dive into the 
details later, for now, just save it to a new file called `somatic_caller.jst`

<details>
<summary>Click here to see the code</summary>

```yaml
# Simple template for calling somatic variants from fastqs for matched 
# tumor/normal pairs
{% set tumor_bam %}{{ tumor.name }}.bam{% endset %}
{% set normal_bam %}{{ normal.nam }}.bam{% endset %}
{% set tumor_rg %}@RG\\tID:{{ tumor.name }}\\tSM:{{ tumor.name }}{% endset %}
{% set normal_rg %}@RG\\tID:{{ normal.name }}\\tSM:{{ normal.name }}{% endset %}

- name: align_tumor
  cmd: |
     set -ue

     # Align fastqs with bwa, sort with samtools
     bwa mem \
       -t {{ threads }} \
       -R "{{ tumor_rg }}" \
       "{{ bwa_index }}" \
       "{{ tumor.r1fq }}" \
       "{{ tumor.r2fq }}" |\
       samtools sort \
         -O BAM \
         -@ {{ threads }} \
         - \
         -o "{{ tumor.name }}.bam"

      # Generate an index with samtools
      samtools index \
        -@ {{ threads }} \
        "{{ tumor.name }}.bam"


- name: align_normal
  cmd: |
     set -ue

     # Align fastqs with bwa, sort with samtools
     bwa mem \
       -t {{ threads }} \
       -R "{{ normal_rg }}" \
       "{{ bwa_index }}" \
       "{{ normal.r1fq }}" \
       "{{ normal.r2fq }}" |\
       samtools sort \
         -O BAM \
         -@ {{ threads }} \
         - \
         -o "{{ normal.name }}.bam"

      # Generate an index with samtools
      samtools index \
        -@ {{ threads }} \
        "{{ normal.name }}.bam"


- name: call_somatic_variants
  cmd: | 
    set -ue

    gatk Mutect2 \
      --reference "{{ reference_fasta }}" \
      --input "{{ tumor_bam }}" \
      --input "{{ normal_bam }}" \
      --tumor-sample "{{ tumor.name }}" \
      --normal-sample "{{ normal.name }}" \
      --output "{{ vcf_path }}"

```
</details>

* In order to run the template, we'll need to provide values to fill in 
the variables. For this example, were going to add all of the input data
to a single config file in YAML format. Save the code below to a new file
called `inputs.yaml`.

<details>
<summary>Click here to see the code</summary>

TODO: These values should be filled in with `${JETSTREAM_TUTORIAL}/<filename>..` 
but first we need to create a tutorial dataset and upload to some public place

```yaml
tumor:
    name: 
    r1fq:
    r2fq:
normal:
    name:
    r1fq:
    r2fq:
bwa_index:
reference_fasta:
```
</details>

* Now that we have a template and an input config file we can run the pipeline.
In order for Jetstream to save progress and run logs, we need to create a new
project directory. When running pipelines inside of a jetstream project directory
logs will be organized, and progress data will allow you to pause/restart/resume
runs. Use the following command to initialize a new project:

```shell
jetstream init js_somatic_tutorial
```

* Before running the pipeline, we'll verify that the template and config data are
valid with the render feature. This option will fill in all of the variables in
the template with data from our config file, and show us the final commands that 
will be executed when the pipeline runs. Use the following command to render the
template (notice the `-C` is a capital c):

```shell
jetstream render somatic_caller.jst -C inputs.yaml
```

* You should have received an error that there was no value provided for the 
variable `threads`. This was intentionally omitted from our config file to 
demonstrate some of the options for providing config data. Several values can be
given with a config file using `-C/--config-file` options, or single config values
can be given with the `-c/--config` options. Since we don't know what kind of 
computer you're using to run this pipeline, so we'll let you decide the number
of threads to use. Here's how you would supply that value via the command line, 
replace 4 with the number of processor cores you would like to use:

```shell
jetstream render somatic_caller.jst -C inputs.yaml -c int:threads 4
```

* Now the template should render successfully and be printed to the console.
If everything looks good, we're ready to run. Use the following command to 
execute the pipeline.

```shell
jetstream run somatic_caller.jst -C inputs.yaml -c int:threads 4 --project js_somatic_tutorial
```

* The runner will start executing commands in the order descibed by the task
dependencies. If any tasks fail, you'll see messages in the logs, and the 
runner will exit with a non-zero status. If all tasks complete successfully
the runner will exit with a zero.

* Since we used a project, all outputs will be saved in that directory. 
After the runner completes, we can `cd js_somatic_tutorial` and see the results.

* Inside the directory you should see the output data files, along with a 
`jetstream` index directory that contains data generated by the runner. There are 
several subcommands that can be used to inspect projects and get information about
the tasks that were executed. Try the following commands:

```
# See details about the project and current status
jetstream project
# See a brief summary for all tasks 
jetstream tasks
# See detailed information for all tasks
jetstream tasks -v
```

* Thank you for following along.

# Installation

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

## Troubleshooting: Command not found after install

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
