import os
import time
import re
import subprocess
import traceback
import logging
from collections import deque
from threading import Thread

import ulid
from jetstream import utils
from .workflow import save

log = logging.getLogger(__name__)


class ThreadWithReturnValue(Thread):
    def __init__(self, group=None, target=None, name=None, args=(), kwargs=None,
                 *, daemon=None):
        Thread.__init__(self, group, target, name, args, kwargs, daemon=daemon)
        self._return = None


    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args, **self._kwargs)
        else:
            raise RuntimeError('Threads must have a target function')

    def join(self, **kwargs):
        Thread.join(self, **kwargs)
        return self._return


def new_run_id():
    run_id = 'js{}'.format(ulid.new().str)
    return run_id


def cleanse_filename(s):
    # shout-out https://stackoverflow.com/a/1007615/3924113
    s = re.sub(r"[^\w\s]", '', s)
    s = re.sub(r"\s+", '_', s)
    return s


def launch(node_id, node_data):
    """Launch a node

    Nodes represent tasks to be completed. This function launches a
    process for the node, and knows how to handle node directives for
    saving stdout/stderr and piping in stdin.
    """
    log.critical('Launching node process: {}'.format(node_id))

    open_fds = []
    result = {
        'return_code': 1,
        'logs': 'Launcher failed to initialize task!',
    }

    try:
        if 'stdout' in node_data:
            out = open(node_data['stdout'], 'w')
            open_fds.append(out)
        else:
            out = subprocess.PIPE

        if 'stderr' in node_data:
            err = open(node_data['stderr'], 'w')
            open_fds.append(err)
        else:
            err = subprocess.STDOUT

        if 'stdin' in node_data:
            stdin = node_data['stdin'].encode()
        else:
            stdin = None

        p = subprocess.Popen(
            node_data['cmd'],
            stdin=subprocess.PIPE,
            stdout=out,
            stderr=err
        )

        stdout, _ = p.communicate(input=stdin)

        if stdout is not None:
            try:
                stdout = stdout.decode()
            except AttributeError as e:
                stdout = 'Unable to decode stdout: {}'.format(e)

        result['logs'] = stdout
        result['return_code'] = p.returncode

        log.critical('Node process exited: {}'.format(node_id))

    except Exception as e:
        log.exception(e)
        result['logs'] = "Launcher failed:\n{}".format(traceback.format_exc())

    finally:
        for fd in open_fds:
            fd.close()

    return result


def _runner(workflow, run_id, run_path):
    log.critical('Runner initialized: {}'.format(run_path))
    workflow_path = os.path.join(run_path, 'workflow.yaml')
    save(workflow, workflow_path)

    tasks = deque()
    while 1:
        try:
            node = next(workflow)
        except StopIteration:
            log.debug('Workflow raised StopIteration')
            break

        if node is None:
            for node_id, result in _handle_completed(tasks, run_id, run_path):
                workflow.__send__(
                    node_id=node_id,
                    return_code=result['return_code'],
                    logs=result['logs']
                )
                save(workflow, workflow_path)
            time.sleep(1)

        else:
            node_id, node_data = node
            log.critical('Runner sending to launch: {}'.format(node_id))

            thread = ThreadWithReturnValue(
                target=launch,
                args=(node_id, node_data)
            )
            thread.start()
            tasks.append((node_id, thread))
            save(workflow, workflow_path)

    log.critical('Run complete!')


def _handle_completed(tasks, run_id, run_path):
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
                    run_path, node_id.replace('\n', '_') + '.log')

                if os.path.exists(log_path):
                    log.warning('{} already exists'.format(log_path))

                with open(log_path, 'w') as fp:
                    print(res['logs'], file=fp)

                yield (node_id, res)

            except TimeoutError:
                tasks.append(next_task)
                log.critical(
                    'Thread join timeout error: {}'.format(node_id))

        next_task = tasks.popleft()

    return


def run_workflow(workflow):
    parent = os.environ.get('JETSTREAM_RUNPATH', '.jetstream')

    run_id = new_run_id()
    run_path = os.path.join(parent, run_id)

    os.environ['JETSTREAM_RUNID'] = run_id
    os.environ['JETSTREAM_RUNPATH'] = run_path

    log.debug('Making new run {}'.format(run_path))
    os.mkdir(run_path)

    with open(os.path.join(run_path, 'created.yaml'), 'w') as fp:
        record = {run_id: utils.fingerprint()}
        utils.yaml.dump(record, stream=fp)

    _runner(workflow, run_id=run_id, run_path=run_path)

