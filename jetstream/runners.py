import os
import time
import subprocess
import traceback
import logging
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
        else:
            raise RuntimeError('Threads must have a target function')

    def join(self, **kwargs):
        Thread.join(self, **kwargs)
        return self._return


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
    def __init__(self, workflow, project=None):
        self.workflow = workflow

        if project is not None:
            self.workflow.project = project

        self._task_queue = deque()
        self._last_save = time.time()

    def _yield(self):
        pass

    def _launch(self, task_id, task):
        """Launch a task

        This function launches the task command as a subprocess, and knows how
        to handle task directives for saving stdout/stderr and piping in stdin
        data."""
        log.debug('Launching task: {}'.format(task_id))

        fds = setup_io_objs(task_id, task, self.workflow.project.log_path)
        rc = 1
        stdout = None
        stderr = None
        logs = None

        try:
            p = subprocess.Popen(
                task.get('cmd') or 'true',
                stdin=subprocess.PIPE,
                stdout=fds['stdout_fd'],
                stderr=fds['stderr_fd'],
                shell=True
            )

            p.communicate(input=fds['stdin_data'])
            rc = p.returncode
            stdout = fds['stdout_path']
            stderr = fds['stderr_path']

            for fd in fds['fds']:
                fd.close()

        except Exception as e:
            log.exception(e)
            rc = 1
            logs = "Error in launcher!\n{}".format(traceback.format_exc())

        finally:
            return rc, logs, stdout, stderr


    def _handle_completed(self):
        """Cycle through the active tasks queue and return completed tasks. """
        sentinel = object()
        self._task_queue.append(sentinel)

        task = self._task_queue.popleft()
        while task is not sentinel:
            task_id, thread = task

            if thread.is_alive():
                self._task_queue.append(task)
            else:
                try:
                    yield task_id, thread.join(timeout=1)
                except TimeoutError:
                    self._task_queue.append(task)
                    log.critical('Thread timeout error: {}'.format(task_id))

            task = self._task_queue.popleft()

        return

    def start(self):
        log.critical('Initialized:\n{}'.format(self.workflow))
        self._yield()

        while 1:
            try:
                task = next(self.workflow)
            except StopIteration:
                log.debug('Workflow raised StopIteration')
                break

            if task is None:
                for task_id, result in self._handle_completed():
                    if result is None:
                        self.workflow.fail(task_id, logs='Thread crash!')
                    else:
                        rc, logs, stdout, stderr = result

                        if rc != 0:
                            self.workflow.fail(
                                task_id,
                                return_code=rc,
                                stdout=stdout,
                                stderr=stderr,
                                logs=logs)
                        else:
                            self.workflow.complete(
                                task_id,
                                return_code=rc,
                                stdout=stdout,
                                stderr=stderr,
                                logs=logs)

            else:
                task_id, node_data = task

                if node_data['cmd'] is None:
                    log.debug('Autocomplete null command: {}'.format(task_id))

                    self.workflow.complete(task_id)

                else:
                    log.debug('Sending to launch: {}'.format(task_id))

                    thread = ThreadWithReturnValue(
                        target=self._launch,
                        name=task_id,
                        args=(task_id, node_data)
                    )

                    thread.start()

                    self._task_queue.append((task_id, thread))

            self._yield()

        fails = [tid for tid, t in self.workflow.tasks(data=True)
                 if t['status'] == 'failed']

        if fails:
            log.critical('\u2620  Some tasks failed! {}'.format(fails))
            return 1
        else:
            log.critical('\U0001F44D Run complete!')
            return 0

    # def start(self):
    #
    #     try:
    #         return self._start()
    #
    #     except Exception:
    #         for t in self._task_queue:
    #             t.
    #         self._task_queue


def setup_io_objs(task_id, task, out_dir):
    res = {'sdtin_data': None, 'fds': set()}

    # Stdout is task stdout path or a log file named by task id
    if 'stdout' in task:
        filename =  utils.cleanse_filename(task['stdout'])
        stdout_path = os.path.join(out_dir, filename)
        stdout_fd = open(stdout_path, 'w')
        res['fds'].add(stdout_fd)
    else:
        default_filename = utils.cleanse_filename(task_id)
        stdout_path = os.path.join(out_dir, default_filename + '.log')
        stdout_fd = open(stdout_path, 'w')
        res['fds'].add(stdout_fd)

    res['stdout_path'] = stdout_path
    res['stdout_fd'] = stdout_fd

    # Stderr is task stderr path or joined with stdout
    if 'stderr' in task:
        filename = utils.cleanse_filename(task['stderr'])
        stderr_path = os.path.join(out_dir, filename)
        stderr_fd = open(stderr_path, 'w')
        res['fds'].add(stderr_fd)
    else:
        stderr_path = stdout_path
        stderr_fd = subprocess.STDOUT

    res['stderr_path'] = stderr_path
    res['stderr_fd'] = stderr_fd

    # Stdin data is encoded from str(stdin) or None
    if 'stdin' in task:
        stdin_data = str(task['stdin']).encode()
    else:
        stdin_data = None

    res['stdin_data'] = stdin_data

    return res


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
