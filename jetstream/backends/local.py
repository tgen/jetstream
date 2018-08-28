from multiprocessing import cpu_count
from asyncio import BoundedSemaphore
from concurrent.futures import CancelledError
from jetstream import log
from jetstream.backends import Backend


def guess_local_cpus(default=4):
    return cpu_count() or default


class LocalBackend(Backend):
    respects = ('cmd', 'stdin', 'stdout', 'stderr', 'cpus')

    def __init__(self, cpus=None):
        """The LocalBackend executes tasks as processes on the local machine.

        :param cpus: If this is None, the number of available CPUs will be
            guessed.
        """
        n_cpus = cpus or guess_local_cpus()
        self._cpu_sem = BoundedSemaphore(n_cpus)
        log.info('LocalBackend initialized with {} cpus'.format(n_cpus))

    def status(self):
        return 'CPUs: {}'.format(self._cpu_sem)

    async def spawn(self, task):
        log.info('Spawn: {}'.format(task))
        log.debug(self.status())

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

            if stdin:
                stdin_fp = open(stdin, 'r')
            else:
                stdin_fp = None

            if stdout:
                stdout_fp = open(stdout, 'w')
            else:
                stdout_fp = None

            if stderr:
                stderr_fp = open(stderr, 'w')
            else:
                stderr_fp = None

            p = await self.subprocess_run_sh(
                cmd, stdin=stdin_fp, stdout=stdout_fp, stderr=stderr_fp
            )

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
