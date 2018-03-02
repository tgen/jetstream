import logging
import os
import time
from collections import deque
from threading import Thread

from jetstream import utils, profile
from jetstream.projects import exc
from jetstream.projects.launchers import default
from jetstream.projects.runs import new_run_id

log = logging.getLogger(__name__)


# TODO: should project functions walk up the directory tree like git?
# see this https://gist.github.com/zdavkeos/1098474 but also consider
# this:
# [rrichholt@dback-login1:~]$ git pull
# fatal: Not a git repository (or any parent up to mount point /home)
# Stopping at filesystem boundary (GIT_DISCOVERY_ACROSS_FILESYSTEM not set).

class ThreadWithReturnValue(Thread):
    def __init__(self, group=None, target=None, name=None, args=(), kwargs=None,
                 *, daemon=None):
        Thread.__init__(self, group, target, name, args, kwargs, daemon=daemon)
        self._return = None

    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args, **self._kwargs)

    def join(self, **kwargs):
        Thread.join(self, **kwargs)
        return self._return


class Project:
    """Internal representation of a project"""
    def __init__(self, path=None):
        if path is None:
            path = './'

        self.name = os.path.basename(path)

        log.critical('Loading project {}'.format(path))

        if not os.path.exists(path):
            raise exc.NotAProject('Path does not exist: {}'.format(path))

        if not os.path.isdir(path):
            raise exc.NotAProject('Path is not a directory: {}'.format(path))

        self.path = path

        target = os.path.join(self.path, profile['RUN_DATA_DIR'])
        if not os.path.exists(target):
            raise exc.NotAProject('Data dir does not exist {}'.format(target))

        target = os.path.join(self.path, profile['RUN_DATA_DIR'])
        if not os.path.isdir(target):
            raise exc.NotAProject('Data dir is not a dir {}'.format(target))

    def serialize(self):
        return {'name': self.name, 'path': self.path}

    @property
    def data_path(self):
        return os.path.join(self.path, 'project.yaml')

    def data(self):
        data = None
        try:
            data = utils.struct(
                action='load',
                format='yaml',
                path=self.data_path
            )
        except FileNotFoundError as err:
            log.warning(err)
        finally:
            return data

    def runs(self):
        run_data_dir = os.path.join(self.path, profile['RUN_DATA_DIR'])
        return os.listdir(run_data_dir)

    def latest_run(self):
        try:
            latest = sorted(self.runs())[-1]
            return latest
        except IndexError:
            return None

    def new_run(self):
        run_id = new_run_id()
        path = os.path.join(self.path, profile['RUN_DATA_DIR'], run_id)
        os.mkdir(path)

        record = {'created': utils.fingerprint()}
        utils.struct(
            action='dump',
            format='yaml',
            obj=record,
            path=os.path.join(path, 'created')
        )

        log.critical('Created new run record {}'.format(run_id))
        return run_id, path

    def run(self, workflow, launcher=default):
        run_id, run_path = self.new_run()

        os.environ['JETSTREAM_RUN_PATH'] = run_path
        os.environ['JESTREAM_RUN_ID'] = run_id

        workflow_path = os.path.join(run_path, 'workflow')

        tasks = deque()
        while 1:
            try:
                task = next(workflow)
            except StopIteration:
                log.critical('Workflow raised StopIteration')
                break

            if task is None:
                for node_id, result in self._find_completed(tasks):
                    workflow.__send__(node_id, result)
                    workflow.save(workflow_path)
                time.sleep(1)
            else:
                node_id, plugin = task
                thread = ThreadWithReturnValue(
                    target=launcher,
                    args=(plugin,)
                )
                thread.start()
                tasks.append((node_id, thread))

        log.critical('Run complete!')

        # for node in workflow.nodes():
        #     if node['return code']
        # TODO this should return something

    def _find_completed(self, tasks):
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

                    if res is None:
                        # This would be due to an exception occuring
                        # in the launcher function.
                        # TODO Is it necessary to recover from these errors?
                        raise RuntimeError(getattr(thread, 'args'))
                    yield (node_id, res)

                except TimeoutError:
                    tasks.append(next_task)
                    log.critical(
                        'Thread join timeout error: {}'.format(node_id))

            next_task = tasks.popleft()

        return


def init():
    os.mkdir('.jetstream')
    log.critical('Initialized project {}'.format(os.getcwd()))

