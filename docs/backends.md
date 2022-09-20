# Backends

Backends are how the runner is configured to execute tasks on various computing
infrastructures. Backends are designed to be extensible. New  backends need only
to implement a coroutine method "spawn". This coroutine should receive a Task
object, execute the task, mark it pass/failed accordingly, then return the task.
They can optionally include a `coroutines` property. This should be a sequence
of asynchronous callables which will be started when the run begins. An example
of this is added in the Slurm backend, where a coroutine is used to monitor the
slurm accounting database and update any jobs that are  running. Note that you
may want to handle the asyncio.CancelledError in these coroutines otherwise an
error will be reported when the runner exits  early. Also, a `cancel` method can
be added to the backend, which will be called when the runner shuts down early.

Here is an example of what a backend needs. Notice the coroutines instead of
standard function defs.

```python

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
```

## Available Backends

__`local`__ - This is the default backend for new installs of jetstream. This operates
on the host machine, and is typically not recommended for the majority of workflows,
but useful for development and testing.

__`local_docker`__ - Similar to the local backend, but here we utilize docker to run
containerized workflows.

__`local_singularity`__ - Similar to the local backend, but here we utilize singularity
to run containerized workflows. However the current expectation is that the container
image is available under a docker/OCI compliant format.

__`slurm`__ - This backend submits tasks to the slurm job scheduler, the majority
of currently developed workflows (Phoenix, Suncity, etc) are designed to work
of a slurm based cluster.

__`slurm_singularity`__ - This is an augmented backend based on the slurm backend.
Here an intermediary runner is introduced, singularity, which allows for running
containerized workflows within a slurm based cluster. The current expectation is
that the container image is available under a docker/OCI compliant format. Also
supports reading locally stored singularity .sif images.

__`dnanexus`__ - A specialized backend for submitting to DNANexus.
