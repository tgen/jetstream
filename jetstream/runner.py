import asyncio
import signal
from datetime import datetime
from contextlib import contextmanager
from asyncio import BoundedSemaphore, Event
import jetstream
from jetstream import utils, log
from jetstream.backends import LocalBackend


@contextmanager
def sigterm_ignored():
    """Context manager that temporarily catches SIGTERM signals, and raises
    a KeyboardInterrupt instead."""

    def signal_handler(signum, frame):
        """Used to replace the default SIGTERM handler"""
        log.debug('SIGTERM received!')
        raise KeyboardInterrupt

    original_sigint_handler = signal.getsignal(signal.SIGTERM)
    signal.signal(signal.SIGTERM, signal_handler)
    try:
        log.debug('Adding SIGTERM handler')
        yield
    finally:
        log.debug('Removing SIGTERM handler')
        signal.signal(signal.SIGTERM, original_sigint_handler)


class Runner:
    def __init__(self, backend=None, max_forks=499, no_task_ready_delay=None,
                 autosave=30, *args, **kwargs):

        if asyncio._get_running_loop() is not None:
            raise RuntimeError("Cannot start from a running event loop")

        self.loop = asyncio.events.new_event_loop()
        asyncio.events.set_event_loop(self.loop)

        if backend is None:
            backend = LocalBackend

        self.autosave_minimum_interval = autosave
        self.backend = backend(self, *args, **kwargs)
        self.no_task_ready_delay = no_task_ready_delay
        self.run_id = None
        self.workflow = None

        # Attributes used internally that should not be modified
        self._conc_sem = BoundedSemaphore(
            value=max_forks or utils.guess_max_forks(),
            loop=self.loop
        )
        self._errs = False
        self._events = []
        self._futures = []

    async def _autosave(self):
        last_save = datetime.now()
        try:
            while 1:
                # Waits for the next task future to return
                e = Event(loop=self.loop)
                self._events.append(e)
                await asyncio.wait_for(e.wait(), timeout=None)

                # If the workflow has been saved too recently, wait until
                # the minimum interval is reached.
                elapsed = (datetime.now() - last_save).seconds
                delay =  self.autosave_minimum_interval - elapsed
                await asyncio.sleep(delay)

                self._save_workflow()
                last_save = datetime.now()
        except asyncio.CancelledError:
            pass

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
        except Exception as e:
            log.exception(e)
            self._halt()
        finally:
            self._conc_sem.release()
            self._futures.remove(future)

    def _save_workflow(self):
        if self.workflow.save_path:
            self.workflow.save()
        elif self.project:
            self.workflow.save(self.project.workflow_file)
        else:
            self.workflow.save(path=f'{self.run_id}.pickle')

    async def _spawn(self, task):
        await self._conc_sem.acquire()
        self.on_spawn(task)
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
                if not task.directives.get('cmd'):
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
        self.loop.create_task(self._autosave())

    def _start_backend_coroutines(self):
        """Starts up any additional coroutines required by the backend"""
        backend_coroutines = getattr(self.backend, 'coroutines', [])

        if not isinstance(backend_coroutines, (list, tuple)):
            raise ValueError('Backend.coroutines must be a list, or tuple')
        else:
            for c in backend_coroutines:
                self.loop.create_task(c())

    async def _yield(self):
        """Since the workflow is a simple synchronous class, it will just
        return None when there are no tasks ready to be launched. This method
        in allows the runner to "pause" checking the workflow until a task
        future returns, or a timeout (whichever happens first). This greatly
        reduces the cpu load of the runner by not constantly checking the
        workflow for new tasks. """
        e = Event(loop=self.loop)
        self._events.append(e)

        try:
            await asyncio.wait_for(e.wait(), timeout=self.no_task_ready_delay)
        except asyncio.TimeoutError:
            pass

    def close(self):
        if self.loop:
            self.loop.close()

    def on_spawn(self, task):
        log.debug('Runner spawn: {}'.format(task))
        try:
            exec_directive = task.directives.get('exec')

            if exec_directive:
                exec(exec_directive, None, {'runner': self, 'task': task})
                self._workflow_iterator = iter(self.workflow)
        except AttributeError:
            pass

    def preflight(self):
        """Called prior to start"""
        log.info(f'Runner preflight procedure')
        self._save_workflow()

    def shutdown(self):
        """Called after shutdown"""
        log.info(f'Runner shutdown procedure')
        self._save_workflow()

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

    def start(self, workflow, project=None, run_id=None):
        """This method is called to start the runner on a workflow."""
        log.debug('Runner starting...')
        started = datetime.now()
        self.workflow = workflow
        self.project = project
        self.run_id = run_id or jetstream.run_id()
        self._errs = False

        if self.no_task_ready_delay is None:
            self.no_task_ready_delay = 0.1 * len(self.workflow)

        with sigterm_ignored():
            self.preflight()
            log.info(f'Starting run: {self.run_id}')

            try:
                self._start_backend_coroutines()
                self._start_autosave()
                self._main = self.loop.create_task(self._spawn_new_tasks())
                self.loop.run_until_complete(self._main)
            except asyncio.CancelledError:
                self._errs = True
            finally:
                self.shutdown()
                futures = asyncio.Task.all_tasks(self.loop)

                if futures:
                    for task in futures:
                        task.cancel()

                    results = self.loop.run_until_complete(
                        asyncio.gather(
                            *futures,
                            loop=self.loop,
                            return_exceptions=True
                        )
                    )
                else:
                    results = []

                for fut, res in zip(futures, results):
                    if isinstance(res, Exception):
                        self._errs = True

                self.loop.close()

                log.info(f'Total run time: {datetime.now() - started}')
                if self._errs:
                    raise RuntimeError('Runner halted unexpectedly!')
