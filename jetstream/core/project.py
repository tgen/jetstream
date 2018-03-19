import logging
import os
import time
import subprocess
from collections import deque
from threading import Thread

from jetstream import exc, utils
from jetstream.core import run
from jetstream.core.settings import profile
from jetstream.plugins import get_plugin

log = logging.getLogger(__name__)

RUN_DATA_DIR = profile['RUN_DATA_DIR']

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

        target = os.path.join(self.path, RUN_DATA_DIR)
        if not os.path.exists(target):
            raise exc.NotAProject('Data dir does not exist {}'.format(target))

        target = os.path.join(self.path, RUN_DATA_DIR)
        if not os.path.isdir(target):
            raise exc.NotAProject('Data dir is not a dir {}'.format(target))

        self._run_id = ''
        self._run_path = ''

    def serialize(self):
        return {'name': self.name, 'path': self.path}

    @property
    def _data_path(self):
        return os.path.join(self.path, 'project.yaml')

    @property
    def _reference_path(self):
        return os.path.join(self.path, 'reference.yaml')

    def reference(self, *args):
        try:
            data = utils.yaml_load(path=self._reference_path)
        except FileNotFoundError as err:
            log.warning(err)
            return None

        if args:
            if len(args) > 1:
                return {k: v for k, v in data['reference'].items() if k in args}
            else:
                return data['reference'][args[0]]
        else:
            return data['meta']

    def meta(self, *args):
        try:
            data = utils.yaml_load(path=self._data_path)
        except FileNotFoundError as err:
            log.warning(err)
            return None

        if args:
            if len(args) > 1:
                return {k: v for k, v in data['meta'].items() if k in args}
            else:
                return data['meta'][args[0]]
        else:
            return data['meta']

    def data(self, **kwargs):
        try:
            data = utils.yaml_load(path=self._data_path)
        except FileNotFoundError as err:
            log.warning(err)
            return None

        if kwargs:
            return utils.filter_documents(data['data'], kwargs)
        else:
            return data['data']

    def runs(self):
        run_data_dir = os.path.join(self.path, RUN_DATA_DIR)
        return os.listdir(run_data_dir)

    def latest_run(self):
        try:
            latest = sorted(self.runs())[-1]
            return latest
        except IndexError:
            return None

    def new_run(self):
        run_id = run.new_run_id()
        path = os.path.join(self.path, RUN_DATA_DIR, run_id)
        os.mkdir(path)

        record = {'created': utils.fingerprint()}
        utils.yaml_dump(
            obj=record,
            path=os.path.join(path, 'created')
        )

        log.critical('Created new run record {}'.format(run_id))
        return run_id, path

    def load_run(self, run_id=None):
        run_dirs = os.path.join(self.path, '.jetstream')

        if run_id is None:
            latest = self.latest_run()
            return run.load_run(os.path.join(run_dirs, latest))
        else:
            run.load_run(os.path.join(run_dirs, run_id))

    def run(self, workflow):
        self._run_id, self._run_path = self.new_run()

        os.environ['JETSTREAM_RUN_PATH'] = self._run_path
        os.environ['JESTREAM_RUN_ID'] = self._run_id
        workflow_path = os.path.join(self._run_path, 'workflow')

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
                node_id, plugin_id = task
                plugin = get_plugin(plugin_id)

                thread = ThreadWithReturnValue(
                    target=execute_plugin,
                    args=(plugin,)
                )
                thread.start()
                tasks.append((node_id, thread))

        log.critical('Run complete!')

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

                    log_path = os.path.join(self._run_path, node_id) + '.log'
                    print('Eventually will write log to ', log_path)
                    print('For now, here\'s the log:\n', res.logs)
                    yield (node_id, res)

                except TimeoutError:
                    tasks.append(next_task)
                    log.critical(
                        'Thread join timeout error: {}'.format(node_id))

            next_task = tasks.popleft()

        return


def execute_plugin(plugin):
    """Launches plugin by guessing the interpreter to from shebang. This
    currently supports Python2/3 and Bash."""
    log.critical('Starting plugin {}'.format(plugin['plugin_id']))

    # We could add explicit shell options to plugins
    # but this just looks for a shebang
    shebang = plugin['script'].splitlines()[0]

    if not shebang.startswith('#!'):
        log.warning('Unable to parse shebang in {}'.format(plugin['plugin_id']))
        shell_cmd = ['bash']
    elif 'python3' in shebang:
        shell_cmd = ['python3']
    elif 'python' in shebang:
        shell_cmd = ['python']
    else:
        shell_cmd = ['bash']

    p = subprocess.Popen(
        shell_cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=os.environ.copy()
    )

    stdout, _ = p.communicate(input=plugin['script'].encode())

    log.critical('Plugin complete {}'.format(plugin['plugin_id']))
    log.critical('Logs\n{}'.format(stdout.decode()))

    result = Result(
        plugin=plugin,
        return_code=p.returncode,
        logs=stdout.decode()
    )

    return result


def init():
    os.makedirs('.jetstream', exist_ok=True)
    log.critical('Initialized project {}'.format(os.getcwd()))

