# Advanced Workflow Design Tips


This page describes some of the more advanced features in Jinja2 that
can be used to write workflow templates.


## Syntax Highlighting

Sublime text works well for writing templates. There is a syntax highlighter 
package for Ansible that will make it easier to read. First install the Ansible
highlighter:

```
⌘-⇧-P
Package Control: Install Package
(Search for Ansible)
```

Then, configure the editor to use the "Ansible" syntax highlighter for `.jst` 
extensions:

`View -> Syntax -> Open all with current extension as ...`

Ansible is a system management utility that uses YAML/Jinja documents very 
similar to workflow templates. Documentation from Ansible covering the YAML 
template syntax can be found at: 
https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html


## Escape from string formatting hell

String formatting with Python has historically been an ugly process. Since the
templating language is built on Python, it suffers the same problems. String 
interpolation (f-strings) in later versions of Python are a really good 
alternative to %s or format method, but unfortunately are not widely supported 
(Jinja2 does not currently support them). To illustrate the problem, here are
some particularly bad examples building up complex strings with the various 
formatting options:

```
{% set cmd = '%s/bwa mem -R "%s" -t %d %s %s %s 2> %s' % (path, read_group, threads, reference, r1_fastq, r2_fastq, work_dir + '/bwa.err') %}
{% set filename = "temp/" + sample.name + "/" + sample.name + "_" + sample.type + exentension %}
{% set rg = "@RG\tID:{}\tPL:{}\tPU:{}\tSM:{}".format(rgid, rgpl, rgpu, rgsm) %}

```

But, the balanced tag form of `set` is here to save the day. Now, you can 
leverage the template syntax to build strings:

```
{% set cmd %}{{ path }}/bwa mem -R "{{ read_group }}" -t {{ threads }} {{ reference }} {{ r1fastq }} {{ r2fastq }} 2> {{ work_dir }}/bwa.err{% endset %}
{% set filename %}temp/{{ sample.name }}/{{ sample.name }}_{{ sample.type }}{{ extension }}{% endset %}
{% set rg %}@RG\tID:{{ rgid }}\tPL:{{ rgpl }}\tPU:{{ rgpu }}\tSM:{{ rgsm }}{% endset %}
```


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
