"""Runners execute workflows.

This allows for customized task runners to be created in the future.
The BaseRunner manages tasks via a queue of threads, and this will
obviously face some performance limits. But, it has the advantage of
virtually zero overhead or delay for monitoring status changes.

When designing a new runner class, keep in mind:

Task directives are not strictly enforced by the workflow class. There is a
schema in jetstream.tasks.spec, but the workflow only checks that the task has
a unique ID. This allows the Task class to decide which directives are
essential, and how they should affect behavior. Runners rely the api described
in BaseTask.
"""

import time
import signal
import logging
from collections import deque
from threading import Event
from jetstream import utils
from jetstream import LocalTask, SlurmTask

log = logging.getLogger(__name__)


class WorkflowRunner(object):
    task_types = {
        None: LocalTask,
        'local': LocalTask,
        'slurm': SlurmTask
    }

    def __init__(self, workflow, max_concurrency=100000, logging_frequency=30,
                 launch_error_penalty=60):
        self.task_queue = deque()
        self.workflow = workflow
        self.max_concurrency = max_concurrency
        self.logging_frequency = logging_frequency
        self.launch_error_penalty = launch_error_penalty
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
        if time.time() - self._last_update > self.logging_frequency:
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

    def _launch_task(self, task_id, task_directives):
        task_type = task_directives.get('type')
        t = self.task_types[task_type](task_id, task_directives)
        success = t.launch()

        if not success:
            raise Exception('Task failed to launch: {}'.format(t))

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
        max_new_submissions = self.max_concurrency - len(self.task_queue)
        for i in range(max_new_submissions):
            if self._kill.is_set():
                break

            task = next(self.workflow)

            if task is None:
                break
            else:
                task_id, task_directives = task
                
                # TODO REVISIT This seems clunky, and I have a feeling it will
                # cause confusion in the future, but... Tasks that are local
                # (or default type) should have a command. If not, it's
                # probably a symbolic task, only added for coordinating
                # execution dependencies. These tasks can be completed
                # immediately, without the need for launching anything. BUT,
                # other task types with a null command still need to be
                # executed, because they have an "implicit" command. SlurmTasks
                # always launch sbatch for example.

                task_type = task_directives.get('type')
                command = task_directives.get('cmd')
                if task_type in (None, 'local') and command is None:
                    self.workflow.complete(task_id)
                else:
                    self._launch_task(task_id, task_directives)

            self._yield()

    def start(self):
        """ Don't override """
        log.critical('Running workflow...')
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
