import os
import asyncio


class BaseBackend(object):
    """To subclass a backend, just override the "spawn" method with a
    coroutine. max_concurrency can be set to limit the number of jobs
    that a backend will allow to spawn concurrently."""
    runner = None
    semaphore = None
    max_concurrency = -1  # Set a Backend-specific limit on the concurrency

    def start(self, runner):
        self.runner = runner
        self.semaphore = asyncio.BoundedSemaphore(
            value=self.max_concurrency,
            loop=self.runner.loop
        )

    async def coro(self):
        pass

    async def spawn(self, task):
        raise NotImplementedError

    def get_fd_paths(self, task):
        """When working inside project, task outputs will be directed into
        log files inside project.logs_dir. But task.stdout/stderr should
        override this behavior. Backends should use this method to get the
        correct output paths for a task."""
        if self.runner.project:
            if 'stdout' in task.directives:
                stdout = task.directives['stdout']
            else:
                name = task.directives.get('name')
                tid = task.tid

                if name:
                    filename = '{}_{}.out'.format(name, tid)
                else:
                    filename = '{}.out'.format(tid)

                stdout = os.path.join(self.runner.project.logs_dir, filename)

            if 'stderr' in task.directives:
                stderr = task.directives['stderr']
            else:
                stderr = stdout

        else:
            if 'stdout' in task.directives:
                stdout = task.directives['stdout']
            else:
                stdout = None

            if 'stderr' in task.directives:
                stderr = task.directives['stderr']
            else:
                stderr = None

        if 'stdin' in task.directives:
            stdin = task.directives['stdin']
        else:
            stdin = None

        return stdin, stdout, stderr

    async def _spawn(self, task):
        try:
            await self.semaphore.acquire()
            await self.spawn(task)
        finally:
            self.semaphore.release()
