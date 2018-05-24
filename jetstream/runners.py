import os
import time
import signal
import logging
import subprocess
from collections import deque
from threading import Event
from jetstream import utils

log = logging.getLogger(__name__)


class BaseTask(object):
    def __init__(self, task_id, task_data):
        self.task_id = task_id
        self.task_data = task_data
        self.stdout_path = None
        self.stderr_path = None
        self.extras = dict()
        self.returncode = -123

        self.stdin_data = task_data.get('stdin')
        if self.stdin_data is not None:
            self.stdin_data = self.stdin_data.encode()

    def poll(self):
        raise NotImplementedError

    def kill(self):
        raise NotImplementedError

    def wait(self):
        raise NotImplementedError

    def launch(self):
        raise NotImplementedError


class LocalTask(BaseTask):
    """ LocalTasks start a process on the local machine.

    Extras can be used to store information about this task that is only
    applicable to a particular type of task."""
    def __init__(self, task_id, task_data):
        super(LocalTask, self).__init__(task_id, task_data)
        self._launched = False
        self._proc = None
        self.fds = set()
        self._stdout_fd = None
        self._stderr_fd = None
        self._setup_out_paths()
        self._setup_out_fds()

    def _default_stdout(self):
        default_filename = utils.cleanse_filename(self.task_id)
        path = os.path.join('logs', default_filename + '.log')
        return path

    def _setup_out_paths(self):
        self.stdout_path = self.task_data.get('stdout', self._default_stdout())
        self.stderr_path = self.task_data.get('stderr', self.stdout_path)

    def _setup_out_fds(self):
        os.makedirs('logs', exist_ok=True)
        stdout_fd = open(self.stdout_path, 'w')
        self.fds.add(stdout_fd)
        self._stdout_fd = stdout_fd

        if self.stderr_path == self.stdout_path:
            self._stderr_fd = subprocess.STDOUT
        else:
            stderr_fd = open(self.stderr_path, 'w')
            self.fds.add(stderr_fd)
            self._stderr_fd = stderr_fd

    @property
    def proc(self):
        if self._proc is None:
            raise AttributeError('This task has not been launched!')
        else:
            return self._proc

    def poll(self):
        return self.proc.poll()

    def kill(self):
        return self.proc.kill()

    def wait(self):
        self.returncode = self.proc.wait()

        for fd in self.fds:
            fd.close()

        return self.returncode

    def launch(self):
        """Launch this task, returns True if launch successful"""
        self._launched = True

        try:
            p = subprocess.Popen(
                self.task_data.get('cmd') or 'true',
                stdin=subprocess.PIPE,
                stdout=self._stdout_fd,
                stderr=self._stderr_fd,
                shell=True
            )

            self.extras['pid'] = p.pid
            self.extras['args'] = p.args

            if self.stdin_data is not None:
                p.stdin.write(self.stdin_data)
                p.stdin.close()

            self._proc = p
            return True

        except BlockingIOError:
            return False


class BaseRunner(object):
    """Runners execute workflows.

    This allows for customized task runners to be created in the future.
    The BaseRunner manages tasks via a queue of threads, and this will
    obviously face some performance limits. But, it has the advantage of
    virtually zero overhead or delay for monitoring status changes.

    When designing a new runner class, keep in mind:

    Task directives are not formally defined. There is a loose schema in
    jetstream.utils.spec, but this is not currently enforced anywhere
    in the workflow class. This allows the runner class to decide which
    directives are essential, and how to use them. The only universally
    essential task attribute is the "id".

    For example, a runner that uses a task broker back end might allow for
    priority directives, which would be useless to the BaseRunner. Runners
    should not raise errors when tasks include a directive that they don't
    need. They may raise errors when a task lacks an essential directive
    (like the path to direct stdout/err for a local process). But, the door
    is intentionally left open to expand/refine task directives and dream
    up new ways of running them."""

    task_types = {
        'local': LocalTask,
        # 'slurm': SlurmTask
    }

    def __init__(self, workflow,  max_tasks=1000, update_frequency=30):
        self.task_queue = deque()
        self.workflow = workflow
        self._update_frequency = update_frequency
        self._max_tasks = max_tasks
        self._kill = Event()
        self._pause = False
        self._last_update = time.time()
        self._completed = len([s for s in self.workflow.status()
                               if s == 'complete'])

    def kill(self):
        log.critical('Halting run!')
        self._kill.set()

    def _yield(self):
        # TODO Scheduler to manage periodic status reports
        if time.time() - self._last_update > self._update_frequency:
            msg = 'Watching {q_len} tasks. {complete}/{total} tasks completed.'

            log.critical(msg.format(
                q_len=len(self.task_queue),
                complete=self._completed,
                total=len(self.workflow)
            ))

            self._last_update = time.time()

    def _signal_handler(self, sig, frame):
        log.critical('{} received!'.format(signal.Signals(2).name))
        self.kill()

    def _finalize(self):
        log.critical('Finalizing run...')

        try:
            self.workflow.save()
        except ValueError as e:
            log.exception(e)

        fails = [tid for tid, t in self.workflow.tasks(data=True)
                 if t['status'] == 'failed']

        if fails:
            log.critical('\u2620  Some tasks failed! {}'.format(fails))
            return 1
        else:
            log.critical('\U0001F44D Run complete!')
            return 0

    def _handle_completed(self, task):
        """ When a completed task is found, this method is called. """
        task.wait()

        if task.returncode != 0:
            self.workflow.fail(
                task.task_id,
                return_code=task.returncode,
                stdout_path=task.stdout_path,
                stderr_path=task.stderr_path)

        else:
            self.workflow.complete(
                task.task_id,
                return_code=task.returncode,
                stdout_path=task.stdout_path,
                stderr_path=task.stderr_path)

        self._completed += 1

    def _launch_task(self, task_id, task_data):
        task_type = task_data.get('type', 'local')
        t = self.task_types[task_type](task_id, task_data)
        t.launch()

        self.task_queue.append(t)

    def _handle_task_queue(self):
        """ Cycle through the task queue once. """
        sentinel = object()
        self.task_queue.append(sentinel)
        task = self.task_queue.popleft()

        while task is not sentinel:
            if self._kill.is_set():
                task.kill()

            if task.poll() is None:
                self.task_queue.append(task)
            else:
                self._handle_completed(task)
                self._pause = False

            self._yield()
            task = self.task_queue.popleft()

        return

    def _check_for_new_tasks(self):
        max_new_submissions = self._max_tasks - len(self.task_queue)
        for i in range(max_new_submissions):
            if self._kill.is_set():
                break

            task = next(self.workflow)

            if task is None:
                break
            else:
                task_id, task_data = task

                if task_data.get('cmd') is None:
                    self.workflow.complete(task_id)
                else:
                    self._launch_task(task_id, task_data)

            self._yield()

    def start(self):
        """ Don't override """
        log.critical('Initialized!')
        signal.signal(signal.SIGINT, self._signal_handler)

        while 1:
            if self._kill.is_set():
                break

            try:
                self._check_for_new_tasks()
            except StopIteration:
                break

            self._handle_task_queue()

        while len(self.task_queue) > 0:
            self._handle_task_queue()

        return self._finalize()

# import sched
# class SchedRunner(object):
#     def __init__(self, workflow, max_tasks=1000):
#         self.workflow = workflow
#         self.max_tasks = max_tasks
#
#         self._scheduler = sched.scheduler()
#         self.task_queue = deque()
#
#         self._kill = Event()
#         self._pause = False
#         self._launched = 0
#         self._completed = 0
#
#
#     def _launch_task(self, task_id, task_data):
#         task_type = task_data.get('type', 'local')
#         t = BaseRunner.task_types[task_type](task_id, task_data)
#         t.launch()
#
#         self.task_queue.append(t)
#
#     def _check_for_new_tasks(self):
#         n =  self.max_tasks - len(self.task_queue)
#
#         log.critical('Checking for new tasks (max {})...'.format(n))
#
#         for i in range(n):
#             if self._kill.is_set():
#                 break
#
#             task = next(self.workflow)
#
#             if task is None:
#                 break
#             else:
#                 task_id, task_data = task
#
#                 if task_data.get('cmd') is None:
#                     self.workflow.complete(task_id)
#                 else:
#                     self._launch_task(task_id, task_data)
#
#     def log_updates(self):
#         log.critical('Watching {} tasks.'.format(len(self.task_queue)))
#         log.critical('{} tasks launched and {} tasks completed.'.format(
#             self._launched, self._completed))
#         self._scheduler.enter(3, 0, self.log_updates)
#
#     def start(self):
#         self._scheduler.enter(3, 0, self.log_updates)
#         self._scheduler.run(blocking=False)
#
#         while 1:
#             try:
#                 self._check_for_new_tasks()
#                 time.sleep(1)
#             except StopIteration:
#                 break
#
#         print(self._scheduler.queue)



def split_fds(task):
    task_filename = utils.cleanse_filename(task['id'])
    out_path = task_filename + '.out'
    err_path = task_filename + '.err'

    out = open(out_path, 'w')
    err = open(err_path, 'w')

    if 'stdin' in task:
        stdin = str(task['stdin']).encode()
    else:
        stdin = None

    return {'stdout': out, 'stderr': err, 'stdin': stdin}
