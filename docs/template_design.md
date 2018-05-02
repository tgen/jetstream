# Template design notes

## Controlling child templates

Because task ids are globally scoped, task ids in child templates need special care.
Also, connecting dependencies inside child templates can be confusing. Here is my 
reccomendation:

Parent template:

```

- id: something
  cmd: do something

- id: another_thing
  cmd: do something else
  after: something

{% with sub_prefix='child_tasks', after="another_thing", before="continue_here" %}
{% include "child.yaml" %}
{% endwith %}

- id: continue_here
  cmd: do more things


```


Child template:

```
# This task serves as a linker between the parent and child template,
# and it allows for a dependency chain to extend into the child 
# template dynamically. This cooperation between the templates but,
# 

# In linker
- id: {{ sub_prefix }}_begin
  after: {{ after }}
  cmd: echo



- id: {{ sub_prefix }}_task
  cmd: do even more things



# Out linker
- id: {{ sub_prefix }}_end
  before: {{ before }}
  cmd: echo 

```


## Child templates given by variable

Include a child template given by a variable. This could be used to select from
a set of templates that all meet certain criteria. But, setting flow control 
directives in child templates can be confusing. 

```
{% if dna_module is defined %}{% include dna_module %}{% endif %}
```

## Prefixing all task ids in a child template

With will create an variable in an isolated scope 

```
{% with task_prefix='alignment_' %}
{% include dna_module %}
{% endwith %}
```

## Filtering child template content

This will be applied to all content included in the child template
```
{% filter replace('id:', 'ID:') %}
{% include dna_module %}
{% endfilter %}
```


## Accessing variables with partial strings

```
[rrichholt@it5687:~]$ jetstream debug_render -s "{{ hello[kit[1:]] }}" --hello '{"kay": 42, "ey": 24}' --kit hey
24
[rrichholt@it5687:~]$ jetstream debug_render -s "{{ hello[kit[1:]] }}" --hello '{"kay": 42, "ey": 24}' --kit okay
42
```
