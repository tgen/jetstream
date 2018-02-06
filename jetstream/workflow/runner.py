""" Runners operate on workflows to execute plugin components. They make use
of workflow methods .__next__() and .__send__() to get tasks and return results
They also require a strategy, which is the function that will be executed for
each plugin component when it is time to launch that component.
"""
import os
import time
import logging
import shutil
from socket import gethostname
from datetime import datetime
from collections import deque
from threading import Thread

from jetstream import plugins
from jetstream.project import Project
from jetstream.utils import yaml, load_yaml
# TODO add a dump_yaml function to utils

log = logging.getLogger(__name__)

_CURRENT_RUN = None


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


def _save(wf):
    """ This function is called every time a wf.__next__() or wf.__send__() is
    issued, and it saves the current workflow state to the project/run/data"""
    run_data_path = os.path.join('.jetstream', _CURRENT_RUN, 'data')
    lock_file = run_data_path + '.lock'

    data = load_yaml(run_data_path)

    data['active'] = {
        'pid': os.getpid(),
        'hostname': gethostname(),
        'last_update': str(datetime.now()),
        'workflow': wf.serialize()
    }

    with open(lock_file, 'w') as fp:
        yaml.dump(data, fp, default_flow_style=False)

    shutil.move(lock_file, run_data_path)


def _handle(tasks, wf):
    """ Cycle through the active tasks to check for completed threads.

    When a completed thread is found, a process record is created and
    sent back to the workflow as a tuple: (node_id, record).

    """
    sentinel = object()
    tasks.append(sentinel)
    next_task = tasks.popleft()
    while next_task is not sentinel:
        node_id, thread = next_task

        if thread.is_alive():
            tasks.append(next_task)
        else:
            try:
                res = thread.join(timeout=1)
                wf.__send__((node_id, res))
                _save(wf)
            except TimeoutError:
                tasks.append(next_task)
                log.critical('Thread join timeout error: {}'.format(node_id))

        next_task = tasks.popleft()


def _threaded_run(wf, strategy):
    """ Sleepy loop that requests new tasks from the workflow, spawns
    new processes, and handles any active processes. """
    tasks = deque()
    while 1:
        try:
            node = next(wf)
        except StopIteration:
            log.critical('wf raised stop iteration')
            break

        if node is None:
            _handle(tasks, wf)
            time.sleep(1)
        else:
            node_id, node_data = node
            plugin = plugins.get_plugin(node_data['pid'])
            thread = ThreadWithReturnValue(target=strategy, args=(plugin, ))
            thread.start()
            tasks.append((node_id, thread))
            _save(wf)


def run(wf, project, strategy):
    global _CURRENT_RUN

    # Make sure we're working in a project dir
    os.chdir(project)
    p = Project(os.getcwd())

    # Request a new run id
    _CURRENT_RUN = p.new_run()

    # Now go!
    _threaded_run(wf, strategy)


def resume(project, run_id):
    # This will be for resuming a given
    pass

