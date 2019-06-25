# Developers guide


> This guide describes implementation details and is not intended for regular
users.

The jetstream package is built around four main classes: Tasks, Workflows, 
Runners, and Backends.

Tasks - are the fundamental unit of a workflow and contain: the commands that 
  need to be run, the directives that control when it should run, description,
  state, etc. 
 
Workflows - are containers for tasks that can connect them together in a 
  DAG (see `Workflow.graph()`). Workflows are used to save progress for runs by
  pickling the entire object. When a new run is started, it can be merged with 
  the project workflow to only run new tasks, this is accomplished with the 
  `mash` function in the workflows module.
  
Runner - There is a single runner class that is used for running workflows. It
  can be configured through application settings, but is mostly tuned for Slurm
  jobs. It basically just works with the WorkflowGraphIterator class to 
  complete every task in a workflow. It is implemented with `asyncio`, the 
  Python asynchronous library, in order to allow a single threaded process that
  seemlessly manages multiple concurrent tasks.

Backend - This is the way the runner is configured to execute tasks on various
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


Task identity
==============

The Task class defines the objects that are added to the DAG. One goal of this 
project was to build an engine that could identify tasks that have already been 
run on a project, and avoid running them again. This was accomplished by 
assigning an identity to every task, which is computed by hashing some of the 
task directives. Task state is not included in the identity so that a "new" 
task can be compared to a "completed". 

Task identity is controlled by the variable: `jetstream.tasks.IDENTITY` which 
is a sequence of the task directives that should be considered part of the 
identity. Later, it might be useful to extend this list to include things like
a docker container id, conda environment, or virtual machine information. For
now, the only directives included are `exec` and `cmd`.

When workflows are mashed together, we first check if the task already exists 
by the task name. If it does exist, then the identities are compared, if the
identity has changed, the task needs to be replaced with the new version, and 
any descendants of that task also need to be reset.


