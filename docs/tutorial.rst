# How to run pipelines

## Profiles

Some settings can be saved in a profile so that they don't need to be changed
every time you run a command. Jetstream will look for a profile json/yaml file
and load settings on startup. The default path to this file is:
`~/.jetstream_profile`, but it can also be set with the environment variable:
`$JETSTREAM_PROFILE`. Here is an example file:

```yaml
backend: slurm
autosave: 60
```

# How to make a pipeline

# 60 seconds to YAML

YAML is a spec that allows you to store common data structures in plain-text

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

Sequence of Mappings Python:

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

Sequence of Mappings YAML:

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


# Congratulations, you're ready

```python
wf = """
- cmd: say Hello, World!

- cmd: hostname > host.txt
  output: host.txt

- cmd: say The host is $(cat host.txt)
  input: host.txt

"""

with open('workflow.yml', 'w') as fp:
    fp.write(wf)
```


```shell
$ jetstream run workflow.yml
[🌵  jetstream] 2018-08-24 14:02:33,062: Version 0.0.10.dev0
[🌵  workflows] 2018-08-24 14:02:33,068: Building workflow...
[🌵  workflows] 2018-08-24 14:02:33,070: Adding task: <Task(new): e4d24d92cbf4b375f6e23fbcba788424bd4d1b88>
[🌵  workflows] 2018-08-24 14:02:33,070: Adding task: <Task(new): e1f24962aed5479d12aa752b26ba7234b2cf3c75>
[🌵  workflows] 2018-08-24 14:02:33,070: Adding task: <Task(new): a7111215cc7d92c77fd70935e690ad56948248b7>
[🌵  workflows] 2018-08-24 14:02:33,071: Workflow ready: <jetstream.Workflow Counter({'new': 3})>
[🌵     runner] 2018-08-24 14:02:33,071: LocalBackend initialized with 8 cpus
[🌵     runner] 2018-08-24 14:02:33,083: Starting Run ID: js01CNPVXTKKD99V015KK9GCYA3S
[🌵     runner] 2018-08-24 14:02:33,084: Workflow manager started!
[🌵     runner] 2018-08-24 14:02:33,084: Spawn: <Task(pending): e1f24962aed5479d12aa752b26ba7234b2cf3c75>
[🌵     runner] 2018-08-24 14:02:33,088: Spawn: <Task(pending): e4d24d92cbf4b375f6e23fbcba788424bd4d1b88>
[🌵     runner] 2018-08-24 14:02:33,099: Done: <Task(complete): e1f24962aed5479d12aa752b26ba7234b2cf3c75>
[🌵     runner] 2018-08-24 14:02:34,093: Spawn: <Task(pending): a7111215cc7d92c77fd70935e690ad56948248b7>
[🌵     runner] 2018-08-24 14:02:34,683: Done: <Task(complete): e4d24d92cbf4b375f6e23fbcba788424bd4d1b88>
[🌵     runner] 2018-08-24 14:02:38,934: Done: <Task(complete): a7111215cc7d92c77fd70935e690ad56948248b7>
[🌵     runner] 2018-08-24 14:02:39,118: Workflow manager stopped!
[🌵     runner] 2018-08-24 14:02:39,118: Shutting down event loop
[🌵     runner] 2018-08-24 14:02:39,119: Run js01CNPVXTKKD99V015KK9GCYA3S Elapsed: 0:00:06.035958
```
