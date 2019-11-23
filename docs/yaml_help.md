# Understanding YAML

YAML is a spec that allows common data structures to be stored as unambiguous 
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
