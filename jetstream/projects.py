import logging
import os
import time
from datetime import datetime
import filelock
import ulid
import jetstream

log = logging.getLogger(__name__)
INDEX = 'jetstream'
LOGS_DIR = 'logs'
HISTORY_DIR = 'history'
INFO_FILENAME = 'project.yaml'
LOCK_FILENAME = 'pid.lock'
WORKFLOW_FILENAME = 'workflow.pickle'
CONFIG_FILENAME = 'config.yaml'


class ProjectInvalid(Exception):
    """Raised when trying to load a project that is missing something"""
    pass


class PidFileLock(filelock.SoftFileLock):
    """This FileLock subclass adds extra info to the lock file upon acquire"""
    def acquire(self, *args, **kwargs):
        super(PidFileLock, self).acquire(*args, **kwargs)
        with open(self.lock_file, 'w') as fp:
            info = jetstream.utils.Fingerprint()
            jetstream.utils.yaml_dump(info.to_dict(), fp)


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
    def __init__(self, path=None):
        if path is None:
            path = os.getcwd()

        timeout = jetstream.settings['projects']['lock_timeout'].get(int)
        self.path = os.path.abspath(path)
        self.name = os.path.basename(path)
        self.lock = PidFileLock(self.pid_path, timeout=timeout)

    def __repr__(self):
        return f'<Project path={self.path}>'

    def exists(self):
        return os.path.exists(self.info_path)

    def get_info(self):
        return jetstream.utils.load_file(self.info_path)

    def get_config(self):
        if os.path.exists(self.config_path):
            return jetstream.utils.load_file(self.config_path)
        else:
            return {}

    def get_history(self):
        yield self.get_info()
        for f in os.listdir(self.history_dir):
            try:
                path = os.path.join(self.history_dir, f)
                yield jetstream.utils.load_yaml(path)
            except jetstream.utils.yaml.YAMLError:
                log.exception(f'Failed to load history file: {path}')

    def get_workflow(self):
        return jetstream.load_workflow(self.workflow_path)

    def init(self, id=None, config=None):
        """Sets up the directories for this project, will overwrite
        jetstream/project.yaml if this already exists. """
        fingerprint = jetstream.utils.Fingerprint(note='init').to_dict()
        id_formatter = jetstream.settings['projects']['id_format'].get(str)
        fingerprint['id'] = id or jetstream.guid(id_formatter)

        os.makedirs(self.index_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        os.makedirs(self.history_dir, exist_ok=True)

        with open(self.info_path, 'w') as fp:
            jetstream.utils.yaml_dump(fingerprint, fp)

        if config is not None:
            self.save_config(config)

        return self.info_path

    def save_config(self, config):
        with open(self.config_path, 'w') as fp:
            jetstream.utils.yaml_dump(config, fp)
        self.save_history('Updated config', data=config)

    def save_history(self, note=None, data=None):
        tries = jetstream.settings['projects']['max_tries_history'].get(int)
        format = jetstream.settings['projects']['filename_format_history'].get(
            str)
        fingerprint = jetstream.utils.Fingerprint(note=note).to_dict()
        fingerprint['data'] = data

        os.makedirs(self.history_dir, exist_ok=True)

        for _ in range(0, tries):
            name = datetime.utcnow().strftime(format)
            path = os.path.join(self.history_dir, name)
            try:
                with open(path, 'x') as fp:
                    jetstream.utils.yaml_dump(fingerprint, fp)
                return path
            except FileExistsError:
                time.sleep(1)
        else:
            raise FileExistsError(f'Unable to create history file: {path}')

    @property
    def config_path(self):
        return os.path.join(self.index_dir, CONFIG_FILENAME)

    @property
    def index_dir(self):
        return os.path.join(self.path, INDEX)

    @property
    def info_path(self):
        return os.path.join(self.index_dir,INFO_FILENAME)

    @property
    def logs_dir(self):
        return os.path.join(self.index_dir, LOGS_DIR)

    @property
    def history_dir(self):
        return os.path.join(self.index_dir, HISTORY_DIR)

    @property
    def pid_path(self):
        return os.path.join(self.index_dir, LOCK_FILENAME)

    @property
    def workflow_path(self):
        return os.path.join(self.index_dir, WORKFLOW_FILENAME)


def new_project(path=None, config=None):
    """Helper function for creating a projects when they may already exist"""
    project = Project(path=path)

    if not project.exists():
        project.init()

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

    # Load any existing config data already in the project, merge with
    # current config, and save back to project.
    existing_config = project.get_config()
    new_config = jetstream.utils.config_stack(existing_config, config)
    project.save_config(new_config)
    return project
