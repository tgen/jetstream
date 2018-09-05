import traceback
from datetime import datetime
import asyncio
from asyncio import BoundedSemaphore
from jetstream.workflows import save_workflow
from jetstream import utils, log
from jetstream.backends import LocalBackend


class Runner:
    def __init__(self, autosave=3600, backend=None, loop=None,
                 max_concurrency=None, project=None,workflow=None,):
        self.autosave = autosave
        self.backend = backend or LocalBackend()
        self.loop = loop
        self.max_conc = max_concurrency or utils.guess_max_forks()
        self.project = project
        self.workflow = workflow
        self.fp = None
        self._backend = None
        self._conc_sem = None
        self._errs = set()
        self._last_save = datetime.now()
        self._main_future = None
        self._secondary_future = None
        self._task_futures = set()

    def start(self, workflow=None, project=None, debug=False, loop=None):
        self.fp = utils.Fingerprint()
        self.workflow = workflow or self.workflow
        self.project = project or self.project
        self.loop = loop or self.loop

        if self.workflow is None:
            raise ValueError('No workflow has been assigned')

        if asyncio._get_running_loop() is not None:
            raise RuntimeError("Cannot start() from a running event loop")

        if self.loop is None:
            self.loop = asyncio.events.new_event_loop()
            asyncio.events.set_event_loop(self.loop)
            self.loop.set_debug(debug)

        log.info('Starting run: {}'.format(self.fp.id))

        self.backend.start(self)
        self._conc_sem = BoundedSemaphore(self.max_conc, loop=self.loop)
        self._main_future = self.loop.create_task(self.main(self.workflow))
        self._secondary_future = self.loop.create_task(self.backend.coro())

        try:
            self.loop.run_until_complete(self._main_future)
        finally:
            self.shutdown()

            if self._errs:
                for err in self._errs:
                    print(err)
                raise RuntimeError('Runner was halted unexpectedly!')

            log.warning('Run complete: {}'.format(self.fp.id))

    async def main(self, tasks):
        try:
            for task in tasks:
                if task is None:
                    await self.sleep(1)
                    continue
                else:
                    await self.spawn(task)

                await self.sleep(0)

            while self._task_futures:
                await self.sleep(1)

            await self.sleep(0)
        except asyncio.CancelledError:
            log.warning('Runner.main was cancelled!')
        finally:
            log.debug('Runner.main stopped')

    async def sleep(self, delay=0):
        log.verbose('Yield for {}s'.format(delay))

        now = datetime.now()
        last = self._last_save
        elapsed = (now - last).seconds

        if self.autosave and elapsed > self.autosave:
            self._last_save = datetime.now()

            if self.project:
                path = self.project.workflow_file
            else:
                # TODO move this to a config module
                path = 'jetstream_workflow_{}.yaml'.format(self.fp.id)

            save_workflow(self.workflow, path=path)

        await asyncio.sleep(delay)

    async def spawn(self, task, *args):
        log.debug('Runner spawn: {}'.format(task))
        await self._conc_sem.acquire()

        log.debug('Creating async task for: {}'.format(task))
        fut = self.loop.create_task(self.backend.spawn(task))
        fut.add_done_callback(self.handle)
        fut.workflow_task = task
        self._task_futures.add(fut)

    def handle(self, future):
        log.debug('Handle: {}'.format(future))
        self._conc_sem.release()

        try:
            future.result()
        except Exception:
            err = 'Unhandled exception in:\n{}\n{}'.format(
                future._coro, traceback.format_exc())
            self.halt(err)
        finally:
            self._task_futures.remove(future)

    def halt(self, err=None):
        if err is not None:
            self._errs.add(err)

        if self._main_future and not self._main_future.cancelled():
            self._main_future.cancel()

    def shutdown(self):
        log.info('Shutting down runner...')
        if not self._secondary_future.done():
            self._secondary_future.cancel()

        to_cancel = asyncio.Task.all_tasks(self.loop)

        if not to_cancel:
            return

        for task in to_cancel:
            task.cancel()

        self.loop.run_until_complete(
            asyncio.gather(*to_cancel, loop=self.loop, return_exceptions=True)
        )

        log.debug('Wrapped up: {}'.format(to_cancel))

        log.debug('Closing event loop...')
        asyncio.events.set_event_loop(None)
        self.loop.close()

#
# class SpawnedTask:
#     def __init__(self, task, backend, loop):
#         self.task = task
#         self.backend = backend
#         self.loop = loop
#         self.future = None
#
#     def spawn(self):
#         self.future = asyncio.ensure_future(
#             self.backend.spawn(self.task),
#             loop=self.loop
#         )
#
#
# class AsyncRunner(object):
#     def __init__(self, backend=None, project=None, workflow=None,
#                  max_forks=None, autosave=0, default_yield=1):
#         """AsyncRunner executes a workflow using a given backend
#         :param backend: Backend object used to spawn tasks
#         :param max_forks: If this is omitted, an estimate of 25% of the system
#             thread limit is used. Over subscribing the fork limit may cause the
#             system to reject creation of new processes, and tasks will fail to
#             launch.
#         """
#         self.backend = backend or LocalBackend()
#         self.project = project
#         self.workflow = workflow
#         self.backend.runner = self
#         self.autosave = timedelta(seconds=autosave)
#         self.default_yield = default_yield
#         self.fp = None
#         self.run = Event()
#         self.tasks = list()
#         self._start_time = None
#         self.loop = None
#         self._workflow_manager = None
#         self._last_save = datetime.now()
#         self._sem = BoundedSemaphore(max_forks or guess_max_forks())
#
#     def log_status(self):
#         """ Logs a status report """
#         log.info('AsyncRunner event-loop load: {}'.format(
#             len(asyncio.Task.all_tasks())))
#         log.info('Workflow status: {}'.format(self.workflow))
#
#     async def _yield(self, delay=None):
#         """Yield is called when the runner can take a break. This allows other
#         coroutines to continue."""
#         if delay is None:
#             delay = self.default_yield
#
#         log.verbose('Yield for {}s'.format(delay))
#
#         if self.autosave and (datetime.now() - self._last_save) > self.autosave:
#             self._last_save = datetime.now()
#
#             if self.project:
#                 path = self.project.workflow_file
#             else:
#                 # TODO move this to a config module
#                 path = 'jetstream_workflow_{}.yaml'.format(self.fp.id)
#
#             save_workflow(self.workflow, path=path)
#
#         await asyncio.sleep(delay)
#
#     async def workflow_manager(self):
#         """Workflow manager will constantly ask workflow for new tasks
#         and spawn a new coroutine for each task when ready. """
#         log.critical('Workflow manager started!')
#
#         try:
#             for task in self.workflow:
#                 if task is None:
#                     await self._yield()
#                 else:
#                     await self.spawn(task)
#                 await self._yield(0)
#         except Exception as e:
#             log.exception(e)
#             log.info('Exception in workflow manager, halting run!')
#             self.loop.stop()
#         finally:
#             log.info('Workflow manager stopped!')
#             self.run.clear()
#
#     async def spawn(self, task):
#         log.verbose('Registering backend.spawn: {}'.format(task))
#         log.verbose('Runner concurrency semaphore: {}'.format(self._sem))
#
#         if task.directives.get('cmd') is None:
#             return task.complete(0)
#
#         try:
#             await self._sem.acquire()
#             asyncio.ensure_future(self.backend.spawn(task))
#         except Exception as e:
#             log.exception(e)
#             log.info('Exception during task spawn, halting run!')
#             self.loop.stop()
#         finally:
#             self._sem.release()
#
#     def start_backend_coros(self):
#         if hasattr(self.backend, 'start_coros'):
#             coros = self.loop.create_task(self.backend.start_coros())
#         else:
#             coros = self.loop.create_task(asyncio.sleep(0))
#
#         return coros
#
#     def finalize_run(self):
#         fails = [t for t in self.workflow.tasks(objs=True) if
#                  t.status == 'failed']
#
#         if fails:
#             log.info('\u2620  Some tasks failed! {}'.format(fails))
#             rc = 1
#         else:
#             rc = 0
#
#         return rc
#
#     def close(self):
#         if self.loop:
#             self.loop.close()
#
#     def start(self, *, workflow=None, project=None, loop=None):
#         self.run.set()
#         self.fp = utils.Fingerprint()
#
#         if workflow:
#             self.workflow = workflow
#
#         if project:
#             self.project = project
#
#         if not self.workflow:
#             raise ValueError('No workflow has been assigned!')
#
#         if self.project:
#             history_file = os.path.join(project.history_dir, self.fp.id)
#             with open(history_file, 'w') as fp:
#                 utils.yaml_dump(self.fp.serialize(), stream=fp)
#
#         start = datetime.now()
#         os.environ.update(JETSTREAM_RUN_ID=self.fp.id)
#         log.info('Starting Run ID: {}'.format(self.fp.id))
#
#         if loop is None:
#             self.loop = asyncio.get_event_loop()
#         else:
#             self.loop = loop
#
#         try:
#             coros = self.start_backend_coros()
#             self.loop.run_until_complete(self.workflow_manager())
#             coros.cancel()
#         except KeyboardInterrupt:
#             log.critical('Received interrupt: Shutting down!')
#             self.run.clear()
#             for t in asyncio.Task.all_tasks():
#                 t.cancel()
#         finally:
#             rc = self.finalize_run()
#             log.info('Shutting down event loop')
#             self.loop.run_until_complete(self.loop.shutdown_asyncgens())
#
#         elapsed = datetime.now() - start
#         log.info('Run {} Elapsed: {}'.format(self.fp.id, elapsed))
#         return rc
#
