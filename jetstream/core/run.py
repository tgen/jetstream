import os
import ulid
import logging
import time
import subprocess
import traceback
from collections import deque
from threading import Thread
from jetstream import utils

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


def launch(node, env):
    log.critical('Launching node {}'.format(node))

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

        if 'stdin' in node:
            #TODO automatically guess if stdin is a file?
            stdin = node['stdin'].encode()
        else:
            stdin = None

        # TODO Had to add this in order to get stuff that's in PATH
        current_env = os.environ.copy()
        current_env.update(env)

        p = subprocess.Popen(
            node['cmd'],
            stdin=subprocess.PIPE,
            stdout=out,
            stderr=err,
            env=current_env,
        )

        stdout, _ = p.communicate(input=stdin)

        try:
            stdout = stdout.decode()

            # TODO How can we get stdout to write to a file and to logs in
            # real-time as the command is running?

        except AttributeError:
            pass

        result['logs'] = stdout
        result['return_code'] = p.returncode

        log.critical('Node complete {}'.format(node))
    except Exception as e:
        log.exception(e)
        result['logs'] = "Launcher failed:\n{}".format(traceback.format_exc())
    finally:
        for fd in open_fds:
            fd.close()

    return result


def _runner(workflow, env):
    log.critical('Runner initialized: {}'.format(env))
    workflow.save(env['JETSTREAM_WORKFLOWPATH'])

    tasks = deque()
    while 1:
        try:
            node = next(workflow)
        except StopIteration:
            log.debug('Workflow raised StopIteration')
            break

        if node is None:
            for node_id, result in _handle_completed(tasks, env):
                workflow.__send__(
                    node_id=node_id,
                    return_code=result['return_code'],
                    logs=result['logs']
                )
                workflow.save(env['JETSTREAM_WORKFLOWPATH'])
            time.sleep(1)

        else:
            log.critical('Runner sending to launch {}'.format(node))
            node_id, node_data = node

            thread = ThreadWithReturnValue(
                target=launch,
                args=(node_data, env)
            )
            thread.start()
            tasks.append((node_id, thread))
            workflow.save(env['JETSTREAM_WORKFLOWPATH'])

    log.critical('Run complete!')


def _handle_completed(tasks, env):
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
                    env['JETSTREAM_RUNPATH'],
                    node_id + '.log'
                )

                if os.path.exists(log_path):
                    log.warning('{} already exists'.format(log_path))

                with open(log_path, 'w') as fp:
                    fp.write(res['logs'])

                yield (node_id, res)

            except TimeoutError:
                tasks.append(next_task)
                log.critical(
                    'Thread join timeout error: {}'.format(node_id))

        next_task = tasks.popleft()

    return


def run_workflow(workflow, project):
    # Generate the environment variables which will be
    # available to each node command.
    run_id = new_run_id()
    run_path = os.path.join(project.path, '.jetstream', run_id)
    log.debug('Making new run {}'.format(run_path))
    os.makedirs(run_path)

    with open(os.path.join(run_path, 'created'), 'w') as fp:
        fp.write(utils.yaml_dumps(utils.fingerprint()))

    env = {
        'JETSTREAM_RUNID': run_id,
        'JETSTREAM_RUNPATH': run_path,
        'JETSTREAM_WORKFLOWPATH': os.path.join(run_path, 'workflow')
    }

    _runner(workflow, env)



