import os

class BaseBackend(object):
    """To subclass a backend, just override the "spawn" method with a
    coroutine. max_concurrency can be set to limit the number of jobs
    that a backend will allow to spawn concurrently."""
    def __init__(self, runner):
        self.runner = runner

    async def coro(self):
        pass

    async def spawn(self, task):
        raise NotImplementedError

    def get_fd_paths(self, task):
        """When working inside project, task outputs will be directed into
        log files inside project.logs_dir. But task.stdout/stderr should
        override this behavior. Backend subclasses should use this method to
        get the correct output paths for a task."""
        stdin = task.directives.get('stdin')

        if self.runner.project:
            if 'stdout' in task.directives:
                stdout = task.directives['stdout']
            else:
                filename = f'{task.label}.log'
                stdout = os.path.join(self.runner.project.logs_dir, filename)

            if 'stderr' in task.directives:
                stderr = task.directives['stderr']
            else:
                stderr = stdout
        else:
            stdout = task.directives.get('stdout')
            stderr = task.directives.get('stderr')

        return stdin, stdout, stderr
