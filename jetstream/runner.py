import asyncio
import logging
import os
import signal
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
    def __init__(self, backend=None, max_forks=None, throttle=None,
                 autosave_min=None, autosave_max=None, *args, **kwargs):
        self.autosave_min = autosave_min or settings['runner']['autosave_min'].get()
        self.autosave_max = autosave_max or settings['runner']['autosave_max'].get()
        self.throttle = throttle or settings['runner']['throttle'].get()
        self.max_forks = max_forks or settings['runner']['max_forks'].get()
        self._conc_sem = None
        self._errs = False
        self._events = []
        self._futures = []
        self._loop = None
        self._last_save = None
        self._workflow_len = None
        self._previous_directory = None
        self._run_started = None

        # Properties
        self._run_id = None
        self._workflow = None
        self._project = None


        if backend is None:
            backend_name = jetstream.settings['runner']['backend'].get()
            backend = jetstream.utils.dynamic_import(backend_name)
        elif isinstance(backend, str):
            backend = jetstream.settings['']

        self._start_event_loop()
        self.backend = backend(self, *args, **kwargs)

    @property
    def run_id(self):
        return self._run_id

    @property
    def workflow(self):
        return self._workflow

    @property
    def project(self):
        return self._project

    @property
    def loop(self):
        return self._loop

    async def _autosave_coro(self):
        self._last_save = datetime.now()
        try:
            while 1:
                mn = timedelta(seconds=self.autosave_min)
                mx = timedelta(seconds=self.autosave_max)

                # If the workflow has been saved too recently, wait until
                # the minimum interval is reached.
                time_since_last = datetime.now() - self._last_save
                delay_for = mn - time_since_last
                await asyncio.sleep(delay_for.seconds)

                # Otherwise, wait for the next future to return and then save
                # but do not wait longer than autosave_max
                e = Event(loop=self.loop)
                self._events.append(e)
                timeout_at = self._last_save + mx
                timeout_in = (timeout_at - datetime.now()).seconds

                try:
                    await asyncio.wait_for(e.wait(), timeout=timeout_in)
                except asyncio.TimeoutError:
                    pass

                self.save_workflow()
                self._last_save = datetime.now()
        except asyncio.CancelledError:
            pass

    def _cleanup_event_loop(self):
        """If the async event loop has outstanding futures, they must be
        cancelled, and results collected, prior to exiting. Otherwise, lots of
        ugly error messages will be shown to the user. """
        futures = asyncio.Task.all_tasks(self.loop)

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

    def _get_workflow_save_path(self):
        if self.workflow.save_path:
            return self.workflow.save_path
        elif self.project:
            return self.project.workflow_file
        else:
            return f'{self.run_id}.pickle'

    def _halt(self):
        if self._main:
            self._main.cancel()

    def _handle_completed_task_futures(self, future):
        """Callback added to task futures when they are spawned"""
        log.debug('Handle: {}'.format(future))
        try:
            future.result()
            if self._events:
                for e in self._events:
                    e.set()
        except Exception:
            log.exception(f'Unhandled exception in a task future: {future}')
            self._halt()
        finally:
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
        for task in self.workflow:
            if task is None:
                await self._yield()
            else:
                self.on_spawn(task)

                if not task.directives().get('cmd'):
                    # if task is blank or None
                    task.complete()
                else:
                    await self._spawn(task)
                    await asyncio.sleep(0)  # Gives event loop a chance to run

        log.info('Workflow complete!')
        if self._futures:
            await asyncio.wait(self._futures)

    def _start_autosave(self):
        """Starts the autosaver coroutine"""
        self.loop.create_task(self._autosave_coro())

    def _start_backend_coroutines(self):
        """Starts up any additional coroutines required by the backend"""
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

    async def _yield(self):
        """Since the workflow is a simple synchronous class, it will just
        return None when there are no tasks ready to be launched. This method
        allows the runner to "pause" checking the workflow until a task
        future returns, or a timeout (whichever happens first). This greatly
        reduces the cpu load of the runner while idle by not constantly
        checking the workflow for new tasks. """
        delay = self.throttle * (self._workflow_len or 0)
        log.debug(f'Yield for {delay}s or when next future returns')
        e = Event(loop=self.loop)
        self._events.append(e)

        try:
            await asyncio.wait_for(e.wait(), timeout=delay)
        except asyncio.TimeoutError:
            pass

    def on_spawn(self, task):
        """Called prior to launching a task"""
        log.debug('Runner spawn: {}'.format(task))
        try:
            exec_directive = task.directives().get('exec')

            if exec_directive:
                exec(exec_directive, None, {'runner': self, 'task': task})
                self._workflow_iterator = iter(self.workflow)
        except AttributeError:
            log.exception('Exception suppressed')

    def preflight(self):
        """Called prior to start"""
        log.info(f'Saving workflow: {self._get_workflow_save_path()}')
        self.save_workflow()
        log.info(f'Starting run: {self.run_id}')

    def save_workflow(self):
        self.workflow.save(self._get_workflow_save_path())

    def shutdown(self):
        """Called after shutdown"""
        log.info(f'Saving workflow: {self._get_workflow_save_path()}')
        self.save_workflow()

        complete = 0
        failed = 0
        for t in self.workflow.tasks(objs=True):
            if t.is_complete():
                complete += 1
            elif t.is_failed():
                failed += 1

        if failed:
            log.info(f'\u2620  {failed} tasks failed!')
            self._errs = True
        else:
            log.info(f'\U0001f44d {complete} tasks complete!')

        log.info(f'Total run time: {datetime.now() - self._run_started}')

    def start(self, project=None, workflow=None, run_id=None):
        """Called to start the runner on a workflow."""
        self._workflow = workflow
        self._project = project
        self._run_id = run_id or jetstream.run_id()
        self._errs = False
        self._run_started = datetime.now()
        self._previous_directory = os.getcwd()
        self._workflow_len = len(workflow)

        if self.max_forks is None:
            self.max_forks = utils.guess_max_forks()

        self._conc_sem = BoundedSemaphore(value=self.max_forks, loop=self.loop)

        with sigterm_ignored():
            self.preflight()

            if self.project:
                os.chdir(self.project.path)

            try:
                self._start_backend_coroutines()
                self._start_autosave()
                self._main = self.loop.create_task(self._spawn_new_tasks())
                self.loop.run_until_complete(self._main)
            except asyncio.CancelledError:
                self._errs = True
            finally:
                self.shutdown()
                self._cleanup_event_loop()
                os.chdir(self._previous_directory)
