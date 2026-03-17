import os
import random


class BaseBackend(object):
    """To subclass a backend, just override the "spawn" method with a
    coroutine. max_concurrency can be set to limit the number of jobs
    that a backend will allow to spawn concurrently."""
    def __init__(self):
        self.runner = None
        self.coroutines = tuple()

    async def spawn(self, task):
        """The base Backend class cannot be used to run workflows, it is
        only for subclassing to make new backends"""
        raise NotImplementedError

    def cancel(self):
        """Called if the run is cancelled for any reason, this is a good
        spot to clean up any outstanding jobs etc.."""
        pass


class PoisonedBackend(BaseBackend):
    async def spawn(self, something):
        if random.random() > 0.95:
            err = 'Unhandled exceptions in the backend should cause the ' \
                  'runner to halt immediately and revoke all jobs.'
            raise Exception(err)
