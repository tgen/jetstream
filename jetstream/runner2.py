import asyncio
import signal
import logging
import random
import traceback
import jetstream
from datetime import datetime
from contextlib import contextmanager
from asyncio import BoundedSemaphore, Event

log = logging.getLogger(__name__)


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
        log.info('Adding SIGTERM handler')
        yield
    finally:
        log.info('Removing SIGTERM handler')
        signal.signal(signal.SIGTERM, original_sigint_handler)


class Backend:
    """
    Backend.coroutines
    Some backends may need to perform duties in addition to spawning tasks.
    Backends can contain a sequence of coroutines that will also be started
    when the runner starts. When the run finishes, they will be cancelled, and
    so they should handle asyncio.CancelledError.
    """

    def __init__(self, runner):
        self.runner = runner
        self.coroutines = tuple()

    async def spawn(self, task):
        raise NotImplementedError(
            'The base Backend class cannot be used to run workflows')


class DummyBackend(Backend):
    def __init__(self, *args, **kwargs):
        super(DummyBackend, self).__init__(*args, **kwargs)
        self.coroutines = (self.coro,)

    async def coro(self):
        try:
            while 1:
                log.info('DummyBackend checking something....')
                await asyncio.sleep(3)
        except asyncio.CancelledError:
            log.info('ok ok im shutting down!')

    async def spawn(self, something):
        try:
            await asyncio.sleep(random.gauss(10,1) + random.random())
            something.complete()
        except KeyboardInterrupt:
            log.info(f'{something} keyboard interrupt!')
        except asyncio.CancelledError:
            log.info(f'{something} cancelled error!')


class PoisonedBackend(Backend):
    async def spawn(self, something):
        try:
            if random.random() > 0.95:
                raise Exception('Bad code in the backend should cause the runner'
                                'to halt!')
        except KeyboardInterrupt:
            log.info(f'{something} keyboard interrupt!')
        except asyncio.CancelledError:
            log.info(f'{something} cancelled error!')


class Runner:
    def __init__(self, backend, max_forks=499, no_task_ready_delay=60,
                 autosave_minimum_interval=5, *args, **kwargs):
        self.backend = backend(self, *args, **kwargs)
        self.autosave_minimum_interval = autosave_minimum_interval
        self.no_task_ready_delay = no_task_ready_delay
        self.run_id = None
        self.workflow = None

        # Attributes used internally that should not be modified
        self._errs = False
        self._futures = []

        if asyncio._get_running_loop() is not None:
            raise RuntimeError("Cannot start from a running event loop")

        self.loop = asyncio.events.new_event_loop()
        asyncio.events.set_event_loop(self.loop)

        self._conc_sem = BoundedSemaphore(max_forks, loop=self.loop)
        self._events = []

    async def _autosave(self):
        last_save = datetime.now()
        try:
            while 1:
                # Schedule an event that will be set when the next task returns
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
        log.debug('Runner spawn: {}'.format(task))
        await self._conc_sem.acquire()

        try:
            exec_directive = task.directives.get('exec')

            if exec_directive:
                exec(exec_directive, None, {'runner': self, 'task': task})
                self._workflow_iterator = iter(self.workflow)

        except AttributeError:
            pass

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
                await self._spawn(task)
                await asyncio.sleep(0)  # Gives other coroutines a chance to run

        log.info('Workflow complete!')
        if self._futures:
            await asyncio.wait(self._futures)

    def _start_autosave(self, min_interval):
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

    def preflight(self):
        """Called prior to start"""
        log.info(f'Runner preflight procedure')
        self._save_workflow()
        log.info(f'Starting run: {self.run_id}')

    def shutdown(self):
        """Called after shutdown"""
        log.info(f'Runner shutdown procedure')
        self._save_workflow()

    def start(self, workflow, project=None, run_id=None, autosave=True):
        """This method is called to start the runner on a workflow."""
        log.debug('Runner starting...')
        started = datetime.now()
        self.workflow = workflow
        self.project = project
        self.run_id = run_id or jetstream.run_id()
        self._errs = False

        with sigterm_ignored():
            self.preflight()

            try:
                self._start_backend_coroutines()
                self._start_autosave(autosave)
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


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    workflow = jetstream.random_workflow()

    r = Runner(backend=DummyBackend)
    r.start(workflow, run_id='banana', autosave=True)
