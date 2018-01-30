""" Async runners operate on workflows to execute tasks. They make use
of workflow methods .__next__() and .__send__() """
import logging
import shlex
import time
import subprocess
from collections import deque

from jetstream import plugins

log = logging.getLogger(__name__)


def spawn(task, **kwargs):
    """ Spawn a subprocess, this returns a coroutine """

    # Capture stdout/stderr by default. This may cause memory issues
    # if the subprocesses produce too much stream data. but less code for now.
    if 'stdout' not in kwargs:
        kwargs['stdout'] = subprocess.PIPE
    if 'stderr' not in kwargs:
        kwargs['stderr'] = subprocess.PIPE

    node_name, node_data = task

    # Prep command for execution on slurm
    cmd = slurm(node_name)
    cmd_args = shlex.split(cmd)

    proc = subprocess.Popen(cmd_args, **kwargs)
    log.critical('Spawned: {} pid: {}'.format(cmd, proc.pid))

    return proc


def handle(tasks, wf):
    """ Cycle through the active tasks to check for completed processes.

    When a completed process is found, a process record is created and
    sent back to the workflow as a tuple: (task_id, record).

    """
    sentinel = object()
    tasks.append(sentinel)
    next_task = tasks.popleft()
    while next_task is not sentinel:
        task, proc = next_task
        rc = proc.poll()

        if rc is None:
            tasks.append(next_task)
        else:
            stdout, stderr = proc.communicate()

            record = {
                "pid": proc.pid,
                "returncode": proc.returncode,
                "stdout": stdout,
                "stderr": stderr,
            }

            wf.__send__((task[0], record))

        next_task = tasks.popleft()


def wf_walker(wf):
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
            handle(tasks, wf)
            time.sleep(1)
        else:
            proc = spawn(task)
            tasks.append((task, proc))


def run(wf, debug=True):
    """ Entry point for running a workflow """
    # Todo configure logging
    log.critical('Starting walker')
    wf_walker(wf)
    log.critical('wf appears to be done')


def slurm(plugin):
    #cmd = "srun {}".format(plugins[plugin])
    cmd = "{}".format(plugins[plugin])
    return cmd

def init_project(path):
    pass


