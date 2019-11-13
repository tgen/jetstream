import logging
import os
import time
from datetime import datetime
import filelock
import jetstream

log = logging.getLogger(__name__)
INDEX_DIR = 'jetstream'
INDEX_FILENAME = 'project.yaml'
LOGS_DIR = 'logs'
HISTORY_DIR = 'history'
PID_FILENAME = 'pid.lock'
WORKFLOW_FILENAME = 'workflow.pickle'


# May 27th 2019 - Decided that the settings template data should be
# dropped for now, 1) it adds a layer of complexity that's almost too much
# for me to explain, 2) and it creates new difficult problems (like the one
# below). Now, the template context just includes project>pipeline>command.
# This still leaves one complicated behavior: The precedence is configured
# so that pipeline values are higher than project values, you cannot edit
# the config.yaml to override pipeline settings in the pipeline.yaml, they
# must be passed as command args. But, this pattern seems to be the best of
# two bad options.

# TODO include project config data saved in the user app settings?
# It seems logical to store user settings template data in the
# project config.yaml when creating the project. But, this sets up some
# odd behaviors. Since project config would be higher priority that user
# general config. Changes to user general config will not be updated in a
# project when re-running workflows later. How could this be handled?


class PidFileLock(filelock.SoftFileLock):
    """This FileLock subclass adds extra info to the lock file upon acquire"""
    def acquire(self, *args, **kwargs):
        super(PidFileLock, self).acquire(*args, **kwargs)
        with open(self.lock_file, 'w') as fp:
            info = jetstream.utils.Fingerprint(pid=True)
            jetstream.utils.dump_yaml(info.to_dict(), fp)


class ProjectPaths:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.index_dir = os.path.join(self.path, INDEX_DIR)
        self.index_path = os.path.join(self.index_dir, INDEX_FILENAME)
        self.logs_dir = os.path.join(self.index_dir, LOGS_DIR)
        self.history_dir = os.path.join(self.index_dir, HISTORY_DIR)
        self.pid_path = os.path.join(self.index_dir, PID_FILENAME)
        self.workflow_path = os.path.join(self.index_dir, WORKFLOW_FILENAME)

    def exists(self):
        return os.path.exists(self.index_path)


class Project:
    """Load a Jetstream project. A project is any directory with an index dir
    and project.yaml: `<project>/jetstream/`.

    <project>/
        .
        ├── jetstream/
        │   ├── logs/
        │   │   ├── <task_name>.log
        │   │   └──     ...
        │   ├── pid.lock (if locked)
        │   ├── project.yaml
        │   └── workflow.pickle
        │
        ├── <project files>
        └──      ...

    """
    def __init__(self, path=None):
        path = path or os.getcwd()
        timeout = jetstream.settings['projects']['lock_timeout'].get(int)
        self.paths = ProjectPaths(path)
        self.lock = PidFileLock(self.paths.pid_path, timeout=timeout)

        try:
            self.index = jetstream.utils.load_yaml(self.paths.index_path)
            self.info = self.index['__project__']
        except FileNotFoundError as e:
            err = f'Project index not found: {e}'
            raise FileNotFoundError(err) from None
        except KeyError:
            err = f'Project index does not contain "__project__" key.'
            raise FileNotFoundError(err) from None

    def __repr__(self):
        return f'<Project path={self.paths.path}>'

    def add_to_history(self, note=None, data=None):
        tries = jetstream.settings['projects']['history_max_tries'].get(int)
        format = jetstream.settings['projects']['history_filename'].get(str)
        fingerprint = jetstream.utils.Fingerprint(note=note).to_dict()
        fingerprint['data'] = data

        os.makedirs(self.paths.history_dir, exist_ok=True)

        for _ in range(0, tries):
            name = datetime.utcnow().strftime(format)
            path = os.path.join(self.paths.history_dir, name)
            try:
                with open(path, 'x') as fp:
                    jetstream.utils.dump_yaml(fingerprint, fp)
                return path
            except FileExistsError:
                time.sleep(1)
        else:
            raise FileExistsError(f'Unable to create history file: {path}')

    def history_iter(self):
        yield self.index
        for f in os.listdir(self.paths.history_dir):
            try:
                path = os.path.join(self.paths.history_dir, f)
                yield jetstream.utils.load_yaml(path)
            except jetstream.utils.yaml.YAMLError:
                log.exception(f'Failed to load history file: {path}')

    @property
    def is_locked(self):
        return self.lock.is_locked

    def list_history(self):
        return list(self.history_iter())

    def load_workflow(self):
        try:
            return jetstream.load_workflow(self.paths.workflow_path)
        except FileNotFoundError:
            return jetstream.Workflow(path=self.paths.workflow_path)

    def update_index(self, data):
        self.index = jetstream.utils.config_stack(self.index, data)
        with open(self.paths.index_path, 'w') as fp:
            jetstream.utils.dump_yaml(self.index, fp)
        self.add_to_history('Updated index data', data=self.index)


def is_project(path):
    paths = ProjectPaths(path)
    return paths.exists()


def init(path=None, config=None, id=None):
    """Sets up a new project dir. This will overwrite jetstream/project.yaml if
    the project already exists. """
    if path is None:
        path = os.getcwd()

    paths = ProjectPaths(path)
    os.makedirs(paths.index_dir, exist_ok=True)
    os.makedirs(paths.logs_dir, exist_ok=True)
    os.makedirs(paths.history_dir, exist_ok=True)

    fingerprint = jetstream.utils.Fingerprint(note='init', id=id).to_dict()
    config = config or dict()
    config.update(__project__=fingerprint)

    with open(paths.index_path, 'w') as fp:
        jetstream.utils.dump_yaml(config, fp)

    wf = jetstream.Workflow()
    wf.save(paths.workflow_path)
    return Project(path)

