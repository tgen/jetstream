import traceback
from datetime import datetime
import asyncio
from asyncio import BoundedSemaphore
from jetstream.workflows import save_workflow
from jetstream import utils, log
from jetstream.backends import LocalBackend


class Runner:
    def __init__(self, autosave=3600, backend=None, loop=None,
                 max_concurrency=None, project=None, workflow=None,):
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
        self._secondary_future.add_done_callback(self.handle_backend_coro)

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
                    if not task.directives.get('cmd'):
                        task.complete()
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

    async def spawn(self, task):
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

            if hasattr(future, 'workflow_task'):
                if not future.workflow_task.is_done():
                    raise ValueError('Backend failed to update task status')

        except Exception:
            err = 'Unhandled exception in:\n{}\n{}'.format(
                future._coro, traceback.format_exc())
            self.halt(err)
        finally:
            self._task_futures.remove(future)

    def handle_backend_coro(self, future):
        try:
            future.result()
        except asyncio.CancelledError:
            pass
        except Exception:
            err = 'Unhandled exception in:\n{}\n{}'.format(
                future._coro, traceback.format_exc())
            self.halt(err)

    def halt(self, err=None):
        if err is not None:
            self._errs.add(err)

        if self._main_future and not self._main_future.cancelled():
            log.info('Halting run!')
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
