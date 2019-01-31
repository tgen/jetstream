import logging
import asyncio
import jetstream
from multiprocessing import cpu_count
from asyncio import BoundedSemaphore, create_subprocess_shell, CancelledError

log = logging.getLogger(__name__)


class LocalBackend(jetstream.backends.BaseBackend):
    def __init__(self, cpus=None, blocking_io_penalty=None):
        """The LocalBackend executes tasks as processes on the local machine.

        This contains a semaphore that limits tasks by the number of cpus
        that they require. It requires that self.runner be set to get the
        event loop, so it's not instantiated until preflight.

        :param cpus: If this is None, the number of available CPUs will be
            guessed. This cannot be changed after starting the backend.
        :param blocking_io_penalty: Delay (in seconds) when a BlockingIOError
            prevents a new process from spawning.
        :param max_concurrency: Max concurrency limit
        """
        super(LocalBackend, self).__init__()
        self._cpu_sem = None
        self.cpus = cpus or jetstream.settings['backends']['local']['cpus'] or self.guess_local_cpus()
        self.blocking_io_penaty = blocking_io_penalty or jetstream.settings['backends']['local']['blocking_io_penalty'].get(int)

    def preflight(self):
        self._cpu_sem = BoundedSemaphore(self.cpus, loop=self.runner.loop)
        log.info(f'LocalBackend initialized with {self.cpus} cpus')

    async def spawn(self, task):
        log.debug('Spawn: {}'.format(task))

        if 'cmd' not in task.directives():
            return task.complete(0)

        cmd = task.directives()['cmd']
        cpus = task.directives().get('cpus', 0)
        cpus_reserved = 0
        open_fps = list()

        if cpus > self.cpus:
            raise RuntimeError('Task cpus greater than available cpus')

        try:
            for i in range(task.directives().get('cpus', 0)):
                await self._cpu_sem.acquire()
                cpus_reserved += 1

            stdin, stdout, stderr = self.get_fd_paths(task)

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

            p = await self.subprocess_sh(
                cmd,
                stdin=stdin_fp,
                stdout=stdout_fp,
                stderr=stderr_fp
            )

            log.info(f'LocalBackend spawned({p.pid}): {task.tid}')
            rc = await p.wait()

            if rc != 0:
                log.info(f'Failed: {task.tid}')
                return task.fail(p.returncode)
            else:
                log.info(f'Complete: {task.tid}')
                return task.complete(p.returncode)
        except CancelledError:
            task.state['err'] = 'Runner cancelled Backend.spawn()'
            return task.fail(-15)
        finally:
            for fp in open_fps:
                fp.close()

            for i in range(cpus_reserved):
                self._cpu_sem.release()

    async def subprocess_sh(
            self, args, *, stdin=None, stdout=None, stderr=None,
            cwd=None, encoding=None, errors=None, env=None,
            loop=None, executable='/bin/bash'):
        """Asynchronous version of subprocess.run

        This will always use a shell to launch the subprocess, and it prefers
        /bin/bash (can be changed via arguments)"""
        log.debug('subprocess_run_sh: {}'.format(args))

        while 1:
            try:
                p = await create_subprocess_shell(
                    args,
                    stdin=stdin,
                    stdout=stdout,
                    stderr=stderr,
                    cwd=cwd,
                    encoding=encoding,
                    errors=errors,
                    env=env,
                    loop=loop,
                    executable=executable
                )
                break
            except BlockingIOError as e:
                log.warning('System refusing new processes: {}'.format(e))
                await asyncio.sleep(self.blocking_io_penaty)

        return p

    def guess_local_cpus(self, default=1):
        return cpu_count() or default

