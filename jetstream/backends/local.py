import asyncio
from subprocess import CompletedProcess
from asyncio import BoundedSemaphore, create_subprocess_shell, CancelledError
from multiprocessing import cpu_count
from jetstream import log
from jetstream.backends import BaseBackend


def guess_local_cpus(default=1):
    return cpu_count() or default


class LocalBackend(BaseBackend):
    respects = ('cmd', 'stdin', 'stdout', 'stderr', 'cpus')

    def __init__(self, cpus=None, blocking_io_penalty=10, max_concurrency=100):
        """The LocalBackend executes tasks as processes on the local machine.

        :param cpus: If this is None, the number of available CPUs will be
            guessed.
        """
        self.cpus = cpus or guess_local_cpus()
        self.max_concurrency = max_concurrency
        self.blocking_io_penaty = blocking_io_penalty

    def start(self, runner):
        super(LocalBackend, self).start(runner)
        self._cpu_sem = BoundedSemaphore(self.cpus, loop=self.runner.loop)
        log.info('LocalBackend initialized with {} cpus'.format(self.cpus))

    async def create_subprocess_shell(self, cmd, **kwargs):
        while 1:
            try:
                return await create_subprocess_shell(cmd, **kwargs)
            except BlockingIOError as e:
                log.warning('System refusing new processes: {}'.format(e))
                log.warning('Retry in 10 seconds.')
                await asyncio.sleep(self.blocking_io_penaty)

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

        return CompletedProcess(args, p.returncode, stdout, stderr)

    def status(self):
        return 'CPUs: {}'.format(self._cpu_sem)

    async def spawn(self, task):
        log.info('Spawn: {}'.format(task))
        cpus_reserved = 0

        try:
            if 'cmd' not in task.directives:
                task.complete(0)
                return

            cmd = task.directives['cmd']

            for i in range(task.directives.get('cpus', 0)):
                await self._cpu_sem.acquire()
                cpus_reserved += 1

            stdin, stdout, stderr = self.get_fd_paths(task)

            open_fps = list()

            if stdin:
                stdin_fp = open(stdin, 'r')
                open_fps.append(stdin_fp)
            else:
                stdin_fp = None

            if stdout:
                stdout_fp = open(stdout, 'w')
                open_fps.append(stdout_fp)
            else:
                stdout_fp = None

            if stderr:
                stderr_fp = open(stderr, 'w')
                open_fps.append(stderr_fp)
            else:
                stderr_fp = None

            p = await self.subprocess_run_sh(
                cmd, stdin=stdin_fp, stdout=stdout_fp, stderr=stderr_fp
            )

            for fp in open_fps:
                fp.close()

            if p.returncode != 0:
                return task.fail(p.returncode)
            else:
                return task.complete(p.returncode)

        except CancelledError:
            return task.fail(-15)

        finally:
            log.info('Done: {}'.format(task))

            for i in range(cpus_reserved):
                self._cpu_sem.release()
