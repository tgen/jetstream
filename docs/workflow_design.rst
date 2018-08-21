This is an example of how to develop workflow components that behave
more like the functions you would find in a traditional programming
language.

First lets make a basic macro. It works just like a function,
allowing for input arguments, and "returning" text to the
template:

The first line starts the declaration, "gatk_whatver" is the
name of the macro. "batch" is some iterable containing a set of
intervals where GATK should operate. "sample_name" is the name
of the sample were running. When called, this macro will produce
a GATK task that runs on all the intervals in the batch.

.. code-block::

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



We can call the macro like this:

.. code-block::

    {{ gatk_whatever(batch=[1,2,3], sample_name='foo') }}


Next is an example of a nested macro, it is a macro that calls
another macro. This macro can be used to group one list into
smaller lists and then call some other function for each batch.
arguments can be passed to the nested function, fn, by including
positional args in a list, or key-word args in a dict.

.. code-block::

{% macro for_batch_n(intervals, fn, size=3, args=list(), kwargs=dict() ) %}

{% set batches = intervals | batch(size) | list %}
{% for batch in batches %}
{{ fn(batch, *args, **kwargs) }}
{% endfor %}

{%- endmacro %}

And here we call the batcher function, and tell it to call
"gatk_whatever" for each batch it produces.

.. code-block::

{{ for_batch_n(intervals, gatk_whatever, size=4, args=[sample_name,]) }}
