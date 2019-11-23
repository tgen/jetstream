# Developers guide

> This guide describes implementation details and is not intended for regular
users.


# Contributions and Releases

## Adding features to the project

To add new features to the project, follow this checklist

- [ ] Create feature branch (if you've been working in dev, and just realized
      you need to move to a feature branch see 
      [this post](https://stackoverflow.com/a/1628584/3924113))
- [ ] Make required changes for feature
- [ ] Update README and doc/ items to reflect current code
  - [ ] Update TOC (just paste whole doc into auto-TOC generator)
- [ ] Update changelog with changes since last release
- [ ] Commit and push
- [ ] Open pull request to **develop** branch


## Creating a new release

Only project adminstrators can create a new release. After all features and
bug fixes have been integrated for a specific release, follow this checklist
to create it:

- [ ] Create release branch
- [ ] Merge in any new features for this new release
- [ ] Bump version in two places:
  - [ ] setup.py
  - [ ] jetstream/__init__.py
- [ ] Open a pull request to **master** branch
- [ ] Create new release and tag on Github


# Design details

The jetstream package is built around four main classes: Tasks, Workflows, 
Runners, and Backends.

Tasks - are the fundamental unit of a workflow and contain: the commands that 
  need to be run, the directives that control when it should run, description,
  state, etc. 
 
Workflows - are containers for tasks that can connect them together in a 
  DAG (see `Workflow.graph`). Workflows are used to save progress for runs by
  pickling the entire object. When a new run is started, it can be merged with 
  the project workflow to only run new tasks, this is accomplished with the 
  `mash` function in the workflows module.
  
Runner - There is a single runner class that is used for running workflows. It
  can be configured through application settings, but is mostly tuned for Slurm
  jobs. It basically just works with the WorkflowGraphIterator class to 
  complete every task in a workflow. It is implemented with `asyncio`, the 
  Python asynchronous library, in order to allow a single threaded process that
  seemlessly manages multiple concurrent tasks.

Backend - Backends are how the runner is configured to execute tasks on various
  computing infrastructures. At the time of this writing there is a Local 
  and Slurm backend. But, it's designed to be extensible. New backends need 
  only to implement a coroutine method "spawn". This coroutine should receive
  a Task object, execute the task, mark it pass/failed accordingly, then return
  the task. They can optionally include a `coroutines` property. This should
  be a sequence of asynchronous callables which will be started when the run
  begins. An example of this is added in the Slurm backend, where a coroutine
  is used to monitor the slurm accounting database and update any jobs that 
  are running. Note that you may want to handle the asyncio.CancelledError in 
  these coroutines otherwise an error will be reported when the runner exits 
  early. Also, a `cancel` method can be added to the backend, which will be 
  called when the runner shuts down early (ie keyboard interrupt). 

Here is an example of what a backend needs. Notice the coroutines instead of
standard function defs.

.. code-block:: python

    async def anothercoro():
        print('BAHHHHH')


    class Backend(BaseBackend):
        def __init__(self, failure_rate=0.5, max_concurrency=2):
            self.max_concurrency = max_concurrency
            self.failure_rate = failure_rate
            self.coroutines = [self.scream_every_3_seconds,]

        async def scream_every_3_seconds(self):
            while 1:
                print('AHHHHHHHHHHHH!', self.runner)
                await anothercoro()
                await asyncio.sleep(3)

        async def spawn(self, task):
            await asyncio.sleep(random.randint(1,3))

            if random.random() < self.failure_rate:
                task.fail()
            else:
                task.complete()

            return task

Other modules are intended to help generate these core objects and organize
work into projects. 


## Task identity

The Task class defines the objects that are added to the DAG. One goal of this 
project was to build an engine that could identify tasks that have already been 
run on a project, and avoid running them again. This was accomplished by 
assigning an identity to every task, which is computed by hashing some of the 
task directives. Task state is not included in the identity.

Task identity is controlled by the variable: `jetstream.tasks.IDENTITY` which 
is a list of the task directives that should be considered part of the 
identity. Later, it might be useful to extend this list to include things like
a docker container id, conda environment, or virtual machine information. For
now, the only directives included are `exec` and `cmd`.

When workflows are mashed together, we first check if the task already exists 
by the task name. If it does exist, then the identities are compared, if the
identity has changed, the task needs to be replaced with the new version, and 
any descendants of that task also need to be reset.



## Templates

Templates are a method for generating Task objects (and ultimately a Workflow) 
from a plain-text documents and configuration data. For any given pipeline
we need to read in some user input data, then generate a set of tasks to 
meet the needs of a particular project. My initial idea was to create a set of 
core classes and functions, that could be imported into Python scripts, and 
to aid in building the Workflow object. In practice, this involved the same 
basic pattern for every script: read data from some file, fill in some skeleton
commands with data from the file. In pure-python this meant a lot of ugly 
string formatting lines and text document parsing. Templates were an attempt to
generalize this process and stop the ugly string formatting code. 

Templates leverage Jinja2 to generate a yaml document describing Tasks in the 
Workflow. Each item in the document will instantiate a new Task object. All the
task objects get added to a Workflow. 

Templates allow for a small amount of logic as well, for example we can create 
a task for each item in some list entirely in the template. But, there is no
parser/lexer/etc inside Jetstream for this syntax. The document rendering
is handled entirely by the Jinja2 template engine. After rendering the tasks
are loaded by a relatively simple YAML parsing operation. This templating
on top of the structured document is remarkably powerful when it comes to 
building complex workflows. 

