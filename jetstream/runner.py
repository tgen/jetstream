import os
import traceback
from datetime import datetime
import asyncio
from asyncio import BoundedSemaphore
from jetstream import utils, log, projects
from jetstream.backends import LocalBackend


class Runner:
    def __init__(self, autosave=3600, backend=None, loop=None,
                 max_forks=None, debug=False):
        self.autosave = autosave
        self.backend = backend or LocalBackend()
        self.debug = debug
        self.fingerprint = None
        self.loop = loop
        self.max_forks = max_forks or utils.guess_max_forks()
        self.project = None
        self.workflow = None
        self._delay = None
        self._conc_sem = None
        self._errs = set()
        self._last_save = datetime.now()
        self._main_future = None
        self._secondary_future = None
        self._task_futures = set()

    def _run_preflight(self):
        log.debug('Running preflight checks...')
        if hasattr(self.backend, 'preflight'):
            try:
                getattr(self.backend, 'preflight')()
            except Exception:
                log.critical('Backend failed preflight checks!')
                raise

        if asyncio._get_running_loop() is not None:
            raise RuntimeError("Cannot start from a running event loop")

        if self.loop is None:
            self.loop = asyncio.events.new_event_loop()
            asyncio.events.set_event_loop(self.loop)
            self.loop.set_debug(self.debug)

        # Backends have a semaphore that controls the backend-specific
        # concurrent task limit. So, the call to backend.start() needs to
        # happen after the event loop is added to the runner.
        self.backend.start(self)

        if self.project is not None:
            history_file = os.path.join(
                self.project.history_dir,
                self.fingerprint.id
            )

            if os.path.exists(history_file):
                err = 'The run_id already exists: {}'.format(history_file)
                raise RuntimeError(err)

            with open(history_file, 'w') as fp:
                fp.write(self.fingerprint.to_yaml())

    def _check_for_save(self):
        now = datetime.now()
        last = self._last_save
        elapsed = (now - last).seconds

        if self.autosave and elapsed > self.autosave:
            self._save_workflow()
            self._last_save = datetime.now()

    def _save_workflow(self):
        if self.workflow.save_path:
            self.workflow.save()
        elif self.project:
            self.workflow.save(self.project.workflow_file)
        else:
            path = '{}.pickle'.format(self.fingerprint.id)
            self.workflow.save(path=path)

    async def _yield(self, delay=None):
        if delay is None:
            delay = self._delay
        log.verbose('Yield for {}s'.format(delay))
        self._check_for_save()
        await asyncio.sleep(delay)

    def close(self):
        """Close the event loop. This runner will not be able to start again."""
        if self.loop:
            self.loop.close()

    async def main(self, workflow):
        self._workflow_iterator = iter(workflow)
        try:
            while 1:
                try:
                    task = next(self._workflow_iterator)
                except StopIteration:
                    break

                if task is None:
                    await self._yield()
                    continue
                else:
                    if not task.directives.get('cmd'):
                        task.complete()
                    else:
                        await self.spawn(task)

                await self._yield(0)

            while self._task_futures:
                await self._yield()

            await self._yield(0)
        except asyncio.CancelledError:
            log.warning('Runner.main was cancelled!')
        finally:
            self._save_workflow()
            log.debug('Runner.main stopped')

    async def spawn(self, task):
        log.debug('Runner spawn: {}'.format(task))
        await self._conc_sem.acquire()

        exec_directive = task.directives.get('exec')

        if exec_directive:
            log.debug('Exec directive:\n{}'.format(exec_directive))

            try:
                exec(exec_directive, None, {'runner': self, 'task': task})
            except Exception:
                tb = traceback.format_exc()
                log.critical(
                    'Unhandled exception in exec directive!\n{}'.format(tb)
                )
                task.state['exception'] = tb
                task.fail()
                self.halt()

            self._workflow_iterator = iter(self.workflow)

        log.debug('Creating async task for: {}'.format(task))
        fut = self.loop.create_task(self.backend.spawn(task))
        fut.add_done_callback(self.handle)
        fut.workflow_task = task
        self._task_futures.add(fut)

    def halt(self, err=None):
        if err is not None:
            self._errs.add(err)

        if self._main_future and not self._main_future.cancelled():
            log.info('Halting run!')
            self._main_future.cancel()

    def handle(self, future):
        log.debug('Handle: {}'.format(future))
        self._conc_sem.release()

        try:
            future.result()
            if hasattr(future, 'workflow_task'):
                if not future.workflow_task.is_done():
                    raise RuntimeError(
                        'Backend failed to update task status! A task future '
                        'was completed without the task being marked as done. '
                        'This will eventually allow a task to be automatically '
                        'restarted, but is not supported yet. The backend must '
                        'complete/fail the task before returning.'
                    )
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

    def start(self, workflow, project=None, run_id=None, debug=False,
              check_workflow=True, loop=None):
        log.debug('Configuring runner...')
        started = datetime.now()

        self.workflow = workflow
        self.project = project
        self.fingerprint = utils.Fingerprint(id=run_id)
        self.loop = loop or self.loop
        self.debug = debug or self.debug
        self._run_preflight()

        log.info('Starting run: {}'.format(self.fingerprint.id))
        self._delay = self._delay or (len(self.workflow) * 0.001)
        self._conc_sem = BoundedSemaphore(self.max_forks, loop=self.loop)
        self._main_future = self.loop.create_task(self.main(self.workflow))
        self._secondary_future = self.loop.create_task(self.backend.coro())
        self._secondary_future.add_done_callback(self.handle_backend_coro)

        try:
            self.loop.run_until_complete(self._main_future)
        finally:
            self.shutdown(check_workflow=check_workflow)
            self.fingerprint = None
            self.project = None
            self.workflow = None

            elapsed = datetime.now() - started
            log.critical('Elapsed: {}'.format(elapsed))

            if self._errs:
                for err in self._errs:
                    log.critical(err)
                raise RuntimeError('Runner halted unexpectedly!')

    def shutdown(self, check_workflow=True):
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

        if check_workflow:
            total = len(self.workflow)
            complete = 0
            failed = 0

            for t in self.workflow.tasks(objs=True):
                if t.is_complete():
                    complete += 1
                elif t.is_failed():
                    failed += 1

            if failed:
                msg = f'\u2620  {failed} tasks failed! ' \
                      f'{complete}/{total} complete!'
                raise RuntimeError(msg)
            else:
                log.info(f'\U0001f44d {complete}/{total} tasks complete!')

