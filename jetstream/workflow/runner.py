""" Runners operate on workflows to execute plugin components. They make use
of workflow methods .__next__() and .__send__() to get tasks and return results
"""
import time
import logging
from collections import deque
from threading import Thread

from jetstream import plugins, strategies

log = logging.getLogger(__name__)


class ThreadWithReturnValue(Thread):
    def __init__(self, group=None, target=None, name=None, args=(), kwargs=None, *, daemon=None):
        Thread.__init__(self, group, target, name, args, kwargs, daemon=daemon)
        self._return = None

    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args, **self._kwargs)

    def join(self, **kwargs):
        Thread.join(self, **kwargs)
        return self._return


# def spawn(task, **kwargs):
#     """ Spawn a subprocess, this returns a Popen object """
#
#     # Capture stdout/stderr by default. This may cause memory issues
#     # if the subprocesses produce too much stream data. but less code for now.
#     if 'stdout' not in kwargs:
#         kwargs['stdout'] = subprocess.PIPE
#     if 'stderr' not in kwargs:
#         kwargs['stderr'] = subprocess.PIPE
#
#     node_name, node_data = task
#
#     # Prep command for execution on slurm
#     cmd = slurm(node_name)
#     cmd_args = shlex.split(cmd)
#
#     proc = subprocess.Popen(cmd_args, **kwargs)
#     log.critical('Spawned: {} pid: {}'.format(cmd, proc.pid))
#
#     return proc


def _handle(tasks, wf):
    """ Cycle through the active tasks to check for completed processes.

    When a completed process is found, a process record is created and
    sent back to the workflow as a tuple: (task_id, record).

    """
    sentinel = object()
    tasks.append(sentinel)
    next_task = tasks.popleft()
    while next_task is not sentinel:
        task, thread = next_task

        if thread.is_alive():
            tasks.append(next_task)
        else:
            try:
                res = thread.join(timeout=1)
                wf.__send__((task, res))
            except TimeoutError:
                tasks.append(next_task)
                log.critical('Thread join timeout error: {}'.format(task))

        next_task = tasks.popleft()


def _threaded_run(wf, strategy=strategies.strategy):
    """ Sleepy loop that requests new tasks from the workflow, spawns
    new processes, and handles any active processes. """
    tasks = deque()
    while 1:
        try:
            task = next(wf)
        except StopIteration:
            log.critical('wf raised stop iteration')
            break

        if task is None:
            _handle(tasks, wf)
            time.sleep(1)
        else:
            plugin = plugins.get_plugin(task)
            thread = ThreadWithReturnValue(target=strategy, args=(plugin, ))
            tasks.append((task, thread))


def run(wf, debug=True):
    """ Entry point for running a workflow """
    log.critical('Starting walker')
    _threaded_run(wf)
    log.critical('wf appears to be done')
