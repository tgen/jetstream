import asyncio
import logging
import os
import signal
import sys
from asyncio import BoundedSemaphore, Event
from contextlib import contextmanager
from datetime import datetime, timedelta

import jetstream
from jetstream import utils, settings

log = logging.getLogger(__name__)


@contextmanager
def sigterm_ignored():
    """Context manager that temporarily catches SIGTERM signals and raises
    a KeyboardInterrupt instead."""
    def signal_handler(signum, frame):
        log.debug(f'SIGNAL received: {signum}!')
        raise KeyboardInterrupt

    original_sigint_handler = signal.getsignal(signal.SIGTERM)
    signal.signal(signal.SIGTERM, signal_handler)
    try:
        log.debug('Replacing SIGTERM handler')
        yield
    finally:
        log.debug('Restoring SIGTERM handler')
        signal.signal(signal.SIGTERM, original_sigint_handler)


class Runner:
    def __init__(self, backend_cls=None, backend_params=None,
                 max_concurrency=None, throttle=None, autosave=True,
                 autosave_min=None, autosave_max=None):
        if backend_cls is None:
            self.backend_cls, self.backend_params = jetstream.lookup_backend('local')
        else:
            self.backend_cls = backend_cls
            self.backend_params = backend_params

        self.backend = None
        self.autosave = autosave
        self.autosave_min = autosave_min or settings['runner']['autosave_min'].get()
        self.autosave_max = autosave_max or settings['runner']['autosave_max'].get()
        self.throttle = throttle or settings['runner']['throttle'].get()
        self.max_concurrency = max_concurrency \
                               or settings['runner']['max_concurrency'].get()
        self._conc_sem = None
        self._condition = None
        self._errs = False
        self._futures = []
        self._main = None
        self._previous_directory = None
        self._run_started = None
        self._workflow_iterator = None
        self._workflow_len = None

        # Properties that should not be modified during a run
        self._loop = None
        self._workflow = None
        self._project = None

    @property
    def workflow(self):
        return self._workflow

    @property
    def project(self):
        return self._project

    @property
    def loop(self):
        return self._loop

    # Synchronization methods for events that should wait on tasks to return
    def notify_waiters(self):
        self._event.set()

    async def wait_for_next_task_future(self, timeout=None):
        self._event.clear()
        await asyncio.wait_for(self._event.wait(), timeout=timeout)

    async def _autosave_coro(self):
        log.debug('Autosaver started!')
        last_save = datetime.now()
        try:
            while 1:
                mn = timedelta(seconds=self.autosave_min)
                mx = timedelta(seconds=self.autosave_max)

                # If the workflow has been saved too recently, wait until
                # the minimum interval is reached.
                time_since_last = datetime.now() - last_save
                delay_for = mn - time_since_last
                log.debug(f'Autosaver hold for {delay_for}')
                await asyncio.sleep(delay_for.seconds)

                # Otherwise, wait for the next future to return and then save
                # but do not wait longer than autosave_max
                timeout_at = last_save + mx
                timeout_in = (timeout_at - datetime.now()).seconds
                log.debug(f'Autosaver next save in {timeout_in}s')

                try:
                    await self.wait_for_next_task_future(timeout=timeout_in)
                    log.debug('Autosaver wait cleared by runner')
                except asyncio.TimeoutError:
                    log.debug('Autosaver wait timed out')

                log.debug('Autosaver saving workflow...')
                self.workflow.save()
                last_save = datetime.now()
        except asyncio.CancelledError:
            pass
        finally:
            log.debug('Autosaver stopped!')

    def _cleanup_event_loop(self):
        """If the async event loop has outstanding futures, they must be
        cancelled, and results collected, prior to exiting. Otherwise, lots of
        ugly error messages will be shown to the user. """
        log.debug('Cleanup event loop')
        if sys.version_info < (3, 7):
            futures = asyncio.Task.all_tasks(self.loop)
        else:
            futures = asyncio.all_tasks(self.loop)

        log.debug(f'{len(futures)} outstanding futures to cancel')
        if futures:
            for task in futures:
                task.cancel()

            results = self.loop.run_until_complete(asyncio.gather(
                *futures,
                loop=self.loop,
                return_exceptions=True
            ))

        else:
            results = []

        for fut, res in zip(futures, results):
            if isinstance(res, Exception):
                log.critical(f'Error collected during shutdown: {res}')
                self._errs = True

        if self.loop:
            self.loop.close()

    def _halt(self):
        if self._main:
            self._main.cancel()

    def _handle_completed_task_futures(self, future):
        """Callback added to task futures when they are spawned"""
        log.debug('Handle: {}'.format(future))
        try:
            task = future.result()
            if task.is_failed():
                self._workflow_iterator.graph.skip_descendants(task)
        except Exception:
            log.exception(f'Unhandled exception in a task future: {future}')
            self._halt()
        finally:
            self.notify_waiters()
            self._conc_sem.release()
            self._futures.remove(future)

    async def _spawn(self, task):
        await self._conc_sem.acquire()
        future = self.loop.create_task(self.backend.spawn(task))
        future.add_done_callback(self._handle_completed_task_futures)
        self._futures.append(future)

    async def _spawn_new_tasks(self):
        """This coroutine is where the runner spends most of its time.
        It returns when the workflow raises StopIteration signalling that
        there are no more tasks to run. """
        for task in self._workflow_iterator:
            if task is None:
                await self._yield()
            else:
                self.on_spawn(task)

                if task.is_done():
                    # if there were errors in the exec directive, we do not
                    # attempt to run the cmd
                    await asyncio.sleep(0)
                elif not task.directives.get('cmd'):
                    # if task is blank or None
                    task.complete()
                else:
                    await self._spawn(task)
                    await asyncio.sleep(0)  # Gives event loop a chance to run

        log.info('Run complete!')
        if self._futures:
            await asyncio.wait(self._futures)

    def _start_autosave(self):
        """Starts the autosaver coroutine"""
        self.loop.create_task(self._autosave_coro())

    def _start_backend(self):
        """Starts up any additional coroutines required by the backend"""
        self.backend = self.backend_cls(**self.backend_params)
        self.backend.runner = self
        backend_coroutines = getattr(self.backend, 'coroutines', [])

        if not isinstance(backend_coroutines, (list, tuple)):
            raise ValueError('Backend.coroutines must be a list, or tuple')
        else:
            for c in backend_coroutines:
                self.loop.create_task(c())

    def _start_event_loop(self):
        if asyncio._get_running_loop() is not None:
            raise RuntimeError("Cannot start from a running event loop")

        self._loop = asyncio.events.new_event_loop()
        asyncio.events.set_event_loop(self._loop)
        self._event = Event()

    async def _yield(self):
        """Since the workflow is a simple synchronous class, it will just
        return None when there are no tasks ready to be launched. This method
        allows the runner to "pause" checking the workflow until a task
        future returns, or a timeout (whichever happens first). This greatly
        reduces the cpu load of the runner while idle by not constantly
        checking the workflow for new tasks. """
        delay = self.throttle * (self._workflow_len or 0)
        log.debug(f'Yield for {delay}s or when next future returns')

        try:
            await self.wait_for_next_task_future(timeout=delay)
        except asyncio.TimeoutError:
            pass

    def on_spawn(self, task):
        """Called prior to launching a task"""
        log.debug('Runner spawn: {}'.format(task))

        exec_directive = task.directives.get('exec')

        if exec_directive:
            try:
                exec(exec_directive, None, {'runner': self, 'task': task})
            except Exception as e:
                task.fail(returncode=1)
                task.state['exec error'] = str(e)
                self._workflow_iterator = iter(self.workflow.graph())
                self._workflow_iterator.graph.skip_descendants(task)
            else:
                self._workflow_iterator = iter(self.workflow.graph())

    def preflight(self):
        """Called prior to start"""
        if self.autosave:
            if self.workflow.path:
                log.info(f'Saving workflow: {self.workflow.path}')
                self.workflow.save()
                self._start_autosave()
            else:
                log.warning(
                    'Autosave is enabled, but no path has been set for this '
                    'workflow. Progress will not be saved.'
                )

    def shutdown(self):
        """Called after shutdown"""
        if self.autosave:
            if self.workflow.path:
                log.info(f'Saving workflow: {self.workflow.path}')
                self.workflow.save()

        log.info(f'Total run time: {datetime.now() - self._run_started}')

    def start(self, workflow, project=None):
        """Called to start the runner on a workflow."""
        if project:
            os.chdir(project.paths.path)

        self._project = project
        self._workflow = workflow
        self._workflow_len = len(workflow)
        self._workflow_iterator = iter(self.workflow.graph())
        self._run_started = datetime.now()
        self._previous_directory = os.getcwd()
        self._errs = False
        self._start_event_loop()
        self._start_backend()

        if self.max_concurrency is None:
            self.max_concurrency = utils.guess_max_forks()

        self._conc_sem = BoundedSemaphore(value=self.max_concurrency)

        with sigterm_ignored():
            self.preflight()

            try:
                self._main = self.loop.create_task(self._spawn_new_tasks())
                self.loop.run_until_complete(self._main)
            except asyncio.CancelledError:
                self._errs = True
            finally:
                self.shutdown()
                self._cleanup_event_loop()
                os.chdir(self._previous_directory)
