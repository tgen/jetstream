import os
import subprocess
from datetime import datetime, timedelta
import asyncio
from asyncio import BoundedSemaphore, Event
from jetstream.workflows import save_workflow
from jetstream import utils, log
from jetstream.backends import LocalBackend


class AsyncRunner(object):
    def __init__(self, backend=None, project=None, workflow=None,
                 max_forks=None, autosave=0, default_yield=1):
        """AsyncRunner executes a workflow using a given backend
        :param backend: Backend object used to spawn tasks
        :param max_forks: If this is omitted, an estimate of 25% of the system
            thread limit is used. Over subscribing the fork limit may cause the
            system to reject creation of new processes, and tasks will fail to
            launch.
        """
        self.backend = backend or LocalBackend()
        self.project = project
        self.workflow = workflow
        self.backend.runner = self
        self.autosave = timedelta(seconds=autosave)
        self.default_yield = default_yield
        self.fp = None
        self.run = Event()
        self._start_time = None
        self._loop = None
        self._workflow_manager = None
        self._last_save = datetime.now()
        self._sem = BoundedSemaphore(max_forks or guess_max_forks())

    def log_status(self):
        """ Logs a status report """
        log.info('AsyncRunner event-loop load: {}'.format(
            len(asyncio.Task.all_tasks())))
        log.info('Workflow status: {}'.format(self.workflow))

    async def _yield(self, delay=None):
        if delay is None:
            delay = self.default_yield

        log.verbose('Yield for {}s'.format(delay))

        if self.autosave and (datetime.now() - self._last_save) > self.autosave:
            self._last_save = datetime.now()

            if self.project:
                path = self.project.workflow_file
            else:
                # TODO move this to a config module
                path = 'jetstream_workflow_{}.yaml'.format(self.fp.id)

            save_workflow(self.workflow, path=path)

        await asyncio.sleep(delay)

    async def workflow_manager(self):
        """Workflow manager will constantly ask workflow for new tasks
        and spawn a new coroutine for each task when ready. """
        log.critical('Workflow manager started!')

        try:
            for task in self.workflow:
                if task is None:
                    await self._yield()
                else:
                    await self.spawn(task)
                await self._yield(0)
        except Exception as e:
            log.exception(e)
            log.info('Exception in workflow manager, halting run!')
            self._loop.stop()
        finally:
            log.info('Workflow manager stopped!')
            self.run.clear()

    async def spawn(self, task):
        log.verbose('Registering backend.spawn: {}'.format(task))
        log.verbose('Runner concurrency semaphore: {}'.format(self._sem))

        if task.directives.get('cmd') is None:
            return task.complete(0)

        try:
            await self._sem.acquire()
            asyncio.ensure_future(self.backend.spawn(task))
        except Exception as e:
            log.exception(e)
            log.info('Exception during task spawn, halting run!')
            self._loop.stop()
        finally:
            self._sem.release()

    def start_backend_coros(self):
        if hasattr(self.backend, 'start_coros'):
            coros = self._loop.create_task(self.backend.start_coros())
        else:
            coros = self._loop.create_task(asyncio.sleep(0))

        return coros

    def finalize_run(self):
        fails = [t for t in self.workflow.tasks(objs=True) if
                 t.status == 'failed']

        if fails:
            log.info('\u2620  Some tasks failed! {}'.format(fails))
            rc = 1
        else:
            rc = 0

        return rc

    def close(self):
        if self._loop:
            self._loop.close()

    def start(self, *, workflow=None, project=None, loop=None):
        self.run.set()
        self.fp = utils.Fingerprint()

        if workflow:
            self.workflow = workflow

        if project:
            self.project = project

        if not self.workflow:
            raise ValueError('No workflow has been given!')

        if self.project:
            history_file = os.path.join(project.history_dir, self.fp.id)
            with open(history_file, 'w') as fp:
                utils.yaml_dump(self.fp.serialize(), stream=fp)

        start = datetime.now()
        os.environ.update(JETSTREAM_RUN_ID=self.fp.id)
        log.info('Starting Run ID: {}'.format(self.fp.id))

        if loop is None:
            self._loop = asyncio.get_event_loop()
        else:
            self._loop = loop

        try:
            manager = self._loop.create_task(self.workflow_manager())
            coros = self.start_backend_coros()
            self._loop.run_until_complete(manager)
            coros.cancel()
        except KeyboardInterrupt:
            log.critical('Received interrupt: Shutting down!')
            self.run.clear()
            for t in asyncio.Task.all_tasks():
                t.cancel()
        finally:
            rc = self.finalize_run()
            log.info('Shutting down event loop')
            self._loop.run_until_complete(self._loop.shutdown_asyncgens())

        elapsed = datetime.now() - start
        log.info('Run {} Elapsed: {}'.format(self.fp.id, elapsed))
        return rc



def guess_max_forks(default=500):
    try:
        res = int(0.25 * int(subprocess.check_output('ulimit -u', shell=True)))
        return res
    except FileNotFoundError as e:
        log.exception(e)
        return default

