import os
import subprocess
import asyncio
import asyncio.subprocess
from collections import defaultdict
from asyncio import create_subprocess_shell
from jetstream import log, settings


class Backend(object):
    """Backends should implement the coroutine method "spawn".

    Backend.spawn will be called by the AsyncRunner when a task is ready
    for execution. It should be asynchronous and return the task_id
    and return code as a tuple: (task_id, returncode)

    Backends should declare a set of task directives that are used when
    executing a task: Backend.respects """
    respects = set()
    runner = None

    def get_output_paths(self, task):
        """When working inside project, task outputs will be directed into
        log files inside project.logs_dir. But task.stdout/stderr will
        override this behavior. Use this method to get the correct ouput paths
        for a task."""
        if self.runner.project:
            if 'stdout' in task.directives:
                stdout = task.directives['stdout']
            else:
                filename = settings['task_out_filename_template']
                params = defaultdict(str, **task.serialize())
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

        return stdout, stderr


    async def create_subprocess_shell(self, cmd, **kwargs):
        while self.runner.run.is_set():
            try:
                return await create_subprocess_shell(cmd, **kwargs)
            except BlockingIOError as e:
                log.info('Unable to start subprocess: {}'.format(e))
                log.info('Retry in 10 seconds.')
                await asyncio.sleep(10)

    async def subprocess_run_sh(
            self, args, *, stdin=None, input=None, stdout=None, stderr=None,
            cwd=None, check=False, encoding=None, errors=None, env=None,
            loop=None, executable='/bin/bash'):
        """Asynchronous version of subprocess.run

        This will always use a shell to launch the subprocess, and it prefers
        /bin/bash (can be changed via arguments)"""
        log.debug('subprocess_run_sh: {}'.format(args))

        if stdin:
            pass
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
