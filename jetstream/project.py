import os
import ulid
import time
import shutil
from collections import deque
from threading import Thread
import logging
from jetstream import utils, plugins

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


class Project:
    """Internal representation of a project used for initiating new and
    resuming old runs."""
    def __init__(self, path=None, run_id=None):
        self._run_id = None
        self._workflow = None

        if path is None:
            path = os.getcwd()

        self.path = path
        assert os.path.exists(self.path) and os.path.isdir(self.path)
        # Jetstream projects must have a .jetstream dir
        assert os.path.exists(self.run_data_dir)

        if run_id is None:
            try:
                self._run_id = self.latest_run()
            except IndexError:
                self.new_run()

    @property
    def run_id(self):
        return self._run_id

    @property
    def project_data_path(self):
        return os.path.join(self.path, 'project.yaml')

    @property
    def project_data(self):
        return utils.load_yaml(self.project_data_path)

    @property
    def run_data_dir(self):
        return os.path.join(self.path, '.jetstream')

    @property
    def run_data_path(self):
        if self._run_id is None:
            raise ValueError('Current run is None')
        return os.path.join(self.run_data_dir, self._run_id)

    @property
    def run_data(self):
        return utils.load_yaml(self.run_data_path)

    def load_run(self, run_id):
        if run_id in self.runs():
            self._run_id = run_id
        else:
            raise ValueError('unknown run id: {}'.format(run_id))

    def runs(self):
        runs = []
        for f in os.listdir(self.run_data_dir):
            if is_run(os.path.join(self.run_data_dir, f)):
                runs.append(f)
        return runs

    def latest_run(self):
        runs = self.runs()
        return sorted(runs)[-1]

    def new_run(self, workflow=None):
        run_id = ulid.new().str
        log.critical('Initializing new run {}'.format(run_id))
        self._run_id = run_id
        self._workflow = workflow
        return run_id

    def save(self):
        if not is_valid_run_id(self._run_id):
            raise ValueError('Run ID, {}, is invalid'.format(self._run_id))

        lock_file = os.path.join(self.run_data_path)

        if self._workflow is not None:
            workflow_data = self._workflow.serialize()
        else:
            workflow_data = None

        run_data = {
            'workflow': workflow_data,
            'last_update': utils.fingerprint()
        }

        with open(lock_file, 'w') as fp:
            utils.yaml.dump(run_data, stream=fp, default_flow_style=False)

        shutil.move(lock_file, self.run_data_path)

    def _handle(self, tasks, workflow):
        """Cycle through the active tasks to check for completed threads.

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

                    if res is None:
                        # This would be due to an exception occuring
                        # in the launcher function.
                        # TODO Is it necessary to recover from these errors?
                        raise RuntimeError(thread.args)

                    workflow.__send__((node_id, res))
                    self.save()

                except TimeoutError:
                    tasks.append(next_task)
                    log.critical(
                        'Thread join timeout error: {}'.format(node_id))

            next_task = tasks.popleft()

    def run(self, workflow, strategy):
        self.new_run(workflow=workflow)

        tasks = deque()
        while 1:
            try:
                node = next(workflow)
            except StopIteration:
                log.critical('Workflow raised StopIteration')
                break

            if node is None:
                self._handle(tasks, workflow)
                time.sleep(1)
            else:
                node_id, plugin = node
                thread = ThreadWithReturnValue(target=strategy, args=(plugin,))
                thread.start()
                tasks.append((node_id, thread))

        log.critical('Run complete!')


def is_run(path):
    """ Returns True if path is a valid run """
    if os.path.isfile(path) and is_valid_run_id(os.path.basename(path)):
        try:
            utils.load_yaml(path)
            return True
        except Exception as err:
            log.exception(err)
            return False
    else:
        return False


def is_valid_run_id(id):
    """ Returns True if id is a valid run id """
    try:
        ulid.from_str(id)
        return True
    except (TypeError, ValueError):
        return False
