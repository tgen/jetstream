import os
import logging
import jetstream

log = logging.getLogger(__name__)


class NotAProject(Exception):
    """Raised when trying to load a project that is missing something"""
    pass


class Project:
    """Load a Jetstream project.

    This class is an internal representation of a project. A project is any
    directory with an index dir: `<project>/jetstream/`. Projects can be
    created automatically by instantiating this class with the "new" argument
    set to True: `jetstream.Project(new=True)`. Creating a project this way will
    set up a best-practices directory structure:

    <project>
        .
        ├── jetstream
        │   ├── config
        │   ├── created
        │   ├── workflow
        │   └── history
        ├── logs
        └── ...

    """
    _INDEX = 'jetstream'
    _CREATED = 'created.yaml'
    _HISTORY = 'history'
    _LOGS = 'logs'
    _PID = 'pid.lock'
    _CONFIG = 'config.yaml'
    _WORKFLOW = 'workflow.pickle'

    def __init__(self, path=None):
        if path is not None:
            self.path = os.path.abspath(path)
        else:
            self.path = os.getcwd()

        if not os.path.exists(self.path):
            raise NotAProject(f'Does not exist: {self.path}')

        if not os.path.exists(self.created_file):
            raise NotAProject(f'File not found: {self.created_file}')

    def __repr__(self):
        return f'<Project path={self.path}>'

    @property
    def index_dir(self):
        return os.path.join(self.path, Project._INDEX)

    @property
    def history_dir(self):
        return os.path.join(self.path, Project._HISTORY)

    @property
    def logs_dir(self):
        return os.path.join(self.path, Project._HISTORY)

    @property
    def created_file(self):
        return os.path.join(self.index_dir, Project._CREATED)

    @property
    def pid_file(self):
        return os.path.join(self.index_dir, Project._PID)

    @property
    def config_file(self):
        return os.path.join(self.index_dir, Project._CONFIG)

    def save_config(self, obj):
        with open(self.config_file, 'w') as fp:
            jetstream.utils.yaml_dump(obj, fp)

    def load_config(self):
        with open(self.config_file, 'r') as fp:
            return jetstream.utils.yaml_loads(fp.read())

    def workflow(self):
        """Load the existing workflow for this project."""
        if os.path.exists(self.workflow):
            return jetstream.workflows.load_workflow(self.workflow)
        else:
            log.warning('No existing workflow found in this project!')
            return None


def new_project(path):
    index_dir = os.path.join(path, Project._INDEX)
    os.makedirs(os.path.join(index_dir, Project._HISTORY), exist_ok=True)
    os.makedirs(os.path.join(index_dir, Project._LOGS), exist_ok=True)

    create_file = os.path.join()
    if os.path.exists(self.created_file):
        log.info('Reinitialized project: {}'.format(self.path))
    else:
        with open(self.created_file, 'w') as fp:
            info = jetstream.utils.Fingerprint()
            jetstream.utils.yaml_dump(info.serialize(), stream=fp)
        log.info('Initialized project: {}'.format(self.path))
