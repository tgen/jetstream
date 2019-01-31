import os
import logging
import filelock
import jetstream

log = logging.getLogger(__name__)
sentinel = object()


class ProjectInvalid(Exception):
    """Raised when trying to load a project that is missing something"""
    pass


class PidFileLock(filelock.SoftFileLock):
    """This FileLock subclass adds extra info to the lock file upon acquire"""
    def acquire(self, *args, **kwargs):
        super(PidFileLock, self).acquire(*args, **kwargs)
        with open(self.lock_file, 'w') as fp:
            info = jetstream.utils.Fingerprint()
            jetstream.utils.yaml_dump(info.serialize(), fp)


class Project:
    """Load a Jetstream project.

    This class is an internal representation of a project. A project is any
    directory with an index dir and project.yaml: `<project>/jetstream/`.

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
    _INDEX = 'jetstream'
    _LOGS_DIR = 'logs'
    _INFO_FILENAME = 'project.yaml'
    _LOCK_FILENAME = 'pid.lock'
    _WORKFLOW_FILENAME = 'workflow.pickle'
    _CONFIG_FILENAME = 'config.yaml'
    def __init__(self, path=None, validate=True):
        self._config = sentinel
        self._info = sentinel
        self._workflow = sentinel

        if path is not None:
            self.path = os.path.abspath(path)
        else:
            self.path = os.getcwd()

        if validate:
            if not os.path.exists(self.path):
                raise ProjectInvalid(f'Does not exist: {self.path}')

            if not os.path.exists(self.info_path):
                raise ProjectInvalid(f'File not found: {self.info_path}')

        self.lock = PidFileLock(
            self.pid_path,
            timeout=jetstream.settings['lock_timeout'].get(int)
        )

    def __repr__(self):
        return f'<Project path={self.path}>'

    @property
    def config_path(self):
        return os.path.join(self.index_dir, Project._CONFIG_FILENAME)

    @property
    def index_dir(self):
        return os.path.join(self.path, Project._INDEX)

    @property
    def info_path(self):
        return os.path.join(self.index_dir, Project._INFO_FILENAME)

    @property
    def logs_dir(self):
        return os.path.join(self.index_dir, Project._LOGS_DIR)

    @property
    def pid_path(self):
        return os.path.join(self.index_dir, Project._LOCK_FILENAME)

    @property
    def workflow_path(self):
        return os.path.join(self.index_dir, Project._WORKFLOW_FILENAME)

    @property
    def config(self):
        if self._config is sentinel:
            try:
                self.load_config()
            except FileNotFoundError:
                self._config = None
        return self._config

    @config.setter
    def config(self, value):
        self._config = value

    @property
    def info(self):
        if self._info is sentinel:
            try:
                self.load_info()
            except FileNotFoundError:
                self._info = None
        return self._info

    @info.setter
    def info(self, value):
        self._info = value

    @property
    def workflow(self):
        """Lazy-loads the workflow file when accessed"""
        if self._workflow is sentinel:
            try:
                self.load_workflow()
            except FileNotFoundError:
                self._workflow = None
        return self._workflow

    @workflow.setter
    def workflow(self, value):
        self._workflow = value

    def init(self, workflow=None, config=None):
        os.makedirs(self.index_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)

        self.info = jetstream.utils.Fingerprint().serialize()
        self.config = config or dict()
        self.workflow = workflow or None
        self.save()
        return self

    def load_config(self):
        """Load the config file for this project"""
        self.config = jetstream.utils.load_file(self.config_path)
        return self.config

    def load_info(self):
        self.info = jetstream.utils.load_file(self.info_path)
        return self.info

    def load_workflow(self):
        """Load the existing workflow for this project. This sets
        project.workflow to the workflow object that was loaded."""
        self.workflow = jetstream.workflows.load_workflow(self.workflow_path)
        return self.workflow

    def save_config(self):
        """Save a config file for this project"""
        log.debug(f'Saving project config: {self.config_path}')
        with open(self.config_path, 'w') as fp:
            jetstream.utils.yaml_dump(self.config, fp)

    def save_info(self):
        log.debug(f'Saving project info: {self.info_path}')
        with open(self.info_path, 'w') as fp:
            jetstream.utils.yaml_dump(self.info, fp)

    def save_workflow(self):
        """Save a workflow in this project"""
        log.debug(f'Saving project workflow: {self.workflow_path}')
        jetstream.workflows.save_workflow(self.workflow, self.workflow_path)

    def save(self):
        self.save_config()
        self.save_info()

        if self.workflow:
            self.save_workflow()


def new_project(path=None, config=None):
    """Helper function for creating a new project"""
    p = Project(path=path, validate=False)
    return p.init(config=config)
