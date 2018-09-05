import os
import subprocess
import asyncio
import asyncio.subprocess
from collections import defaultdict
from asyncio import create_subprocess_shell
from jetstream import log, settings


class BaseBackend(object):
    """To subclass a backend, just override the "spawn" method with a
    coroutine. max_concurrency can be set to limit the number of jobs
    that a backend will allow to spawn concurrently."""
    runner = None
    max_concurrency = -1 # Set a Backend-specific limit on the concurrency

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
                filename = settings['task_out_filename_template']
                params = defaultdict(lambda: 'task', **task.serialize())
                filename = filename.format_map(params)
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


class Backend(object):
    """Backends should implement the coroutine method "spawn".

    Backend.spawn will be called by the AsyncRunner when a task is ready
    for execution. It should be asynchronous and return the task_id
    and return code as a tuple: (task_id, returncode)

    Backends should declare a set of task directives that are used when
    executing a task: Backend.respects """
    respects = set()
    runner = None

    def get_fd_paths(self, task):
        """When working inside project, task outputs will be directed into
        log files inside project.logs_dir. But task.stdout/stderr should
        override this behavior. Backends should use this method to get the
        correct output paths for a task."""
        if self.runner.project:
            if 'stdout' in task.directives:
                stdout = task.directives['stdout']
            else:
                filename = settings['task_out_filename_template']
                params = defaultdict(lambda: 'task', **task.serialize())
                filename = filename.format_map(params)
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


    async def create_subprocess_shell(self, cmd, **kwargs):
        while self.runner.run.is_set():
            try:
                return await create_subprocess_shell(cmd, **kwargs)
            except BlockingIOError as e:
                log.warning('Unable to start subprocess: {}'.format(e))
                log.warning('Retry in 10 seconds.')
                await asyncio.sleep(10)

    async def subprocess_run_sh(
            self, args, *, stdin=None, input=None, stdout=None, stderr=None,
            cwd=None, check=False, encoding=None, errors=None, env=None,
            loop=None, executable='/bin/bash'):
        """Asynchronous version of subprocess.run

        This will always use a shell to launch the subprocess, and it prefers
        /bin/bash (can be changed via arguments)"""
        log.debug('subprocess_run_sh: {}'.format(args))

        if stdin and input:
            raise ValueError('Input and Stdin given to subprocess_run_sh.'
                             'Choose only one.')
        else:
            if input:
                stdin = asyncio.subprocess.PIPE
            else:
                stdin = None

        p = await self.create_subprocess_shell(
            args, stdin=stdin, stdout=stdout, stderr=stderr, cwd=cwd,
            encoding=encoding, errors=errors, env=env,
            loop=loop, executable=executable)

        if input:
            stdout, stderr = await p.communicate(input=input)
        else:
            stdout, stderr = await p.communicate()

        if check and p.returncode != 0:
            raise ChildProcessError(args)

        return subprocess.CompletedProcess(args, p.returncode, stdout, stderr)
