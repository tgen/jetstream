import os
import ulid
import logging
import time
import subprocess
import traceback
from collections import deque
from threading import Thread
from jetstream.core import Project

log = logging.getLogger(__name__)


class ThreadWithReturnValue(Thread):
    def __init__(self, group=None, target=None, name=None, args=(), kwargs=None,
                 *, daemon=None):
        Thread.__init__(self, group, target, name, args, kwargs, daemon=daemon)
        self._return = None


    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args, **self._kwargs)
            # try:
            #     self._return = self._target(*self._args, **self._kwargs)
            # except Exception as e:
            #     self._return = e
        else:
            raise RuntimeError('Threads must have a target function')

    def join(self, **kwargs):
        Thread.join(self, **kwargs)
        return self._return


def new_run_id():
    run_id = 'js{}'.format(ulid.new().str)
    return run_id


def launch(node, run_context):
    log.critical('Starting node {}'.format(node['id']))

    open_fds = []
    result = {
        'return_code': 1,
        'logs': 'Launcher failed!',
    }

    try:
        if 'stdout' in node:
            out = open(node['stdout'], 'w')
            open_fds.append(out)
        else:
            out = subprocess.PIPE

        if 'stderr' in node:
            err = open(node['stderr'], 'w')
            open_fds.append(err)
        else:
            err = subprocess.STDOUT

        p = subprocess.Popen(
            node['cmd'],
            stdin=subprocess.PIPE,
            stdout=out,
            stderr=err,
            env=run_context
        )

        stdout, _ = p.communicate(input=node.get('stdin'))

        try:
            stdout = stdout.decode()
        except AttributeError:
            pass

        result['logs'] = stdout
        result['return_code'] = p.returncode

        log.critical('Node complete {}'.format(node['id']))
    except Exception as e:
        log.exception(e)
        result['logs'] = "Launcher failed:\n{}".format(traceback.format_exc())
    finally:
        for fd in open_fds:
            fd.close()

    return result


def _runner(workflow, run_context):
    workflow.save(run_context['JETSTREAM_WORKFLOWPATH'])

    tasks = deque()
    while 1:
        try:
            task = next(workflow)
        except StopIteration:
            log.critical('Workflow raised StopIteration')
            break

        if task is None:
            for node_id, result in _handle_completed(tasks, run_context):
                workflow.__send__(
                    node_id=node_id,
                    return_code=result['return_code'],
                    logs=result['logs']
                )
                workflow.save(run_context['JETSTREAM_WORKFLOWPATH'])
            time.sleep(1)

        else:
            node_id, node_data = task

            thread = ThreadWithReturnValue(
                target=launch,
                args=(node_data, run_context)
            )
            thread.start()
            tasks.append((node_id, thread))

    log.critical('Run complete!')


def _handle_completed(tasks, run_context):
    """Cycle through the active tasks queue and return completed tasks. """
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

                log_path = os.path.join(
                    run_context['JETSTREAM_RUNPATH'],
                    node_id + '.log'
                )
                print('Eventually will write log to ', log_path)
                print('For now, here\'s the log:\n{}'.format(res['logs']))

                yield (node_id, res)

            except TimeoutError:
                tasks.append(next_task)
                log.critical(
                    'Thread join timeout error: {}'.format(node_id))

        next_task = tasks.popleft()

    return


def run_workflow(workflow):
    # Ensure we are working inside of a valid project
    p = Project()

    # Generate the environment variables which will be
    # available to each node command.
    run_id = new_run_id()
    run_path = os.path.join(p.path, '.jetstream', run_id)
    os.makedirs(run_path)
    env = {
        'JETSTREAM_RUNID': run_id,
        'JETSTREAM_RUNPATH': run_path,
        'JETSTREAM_WORKFLOWPATH': os.path.join(run_path, 'workflow')
    }

    _runner(workflow, env)



