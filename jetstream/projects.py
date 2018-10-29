import os
import jetstream
from jetstream import log, settings

# Project index dir configuration
index_dir = settings['project_index_dir']
history_dir = os.path.join(index_dir, 'history')
logs_dir = settings['project_logs_dir']
results_dir = settings['project_results_dir']
temp_dir = settings['project_temp_dir']
config_dir = settings['project_config_dir']
pid_file = os.path.join(index_dir, 'pid')
config_file = os.path.join(index_dir, 'config')
created_file = os.path.join(index_dir, 'created')
workflow_file = os.path.join(index_dir, 'workflow')


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
        ├── results
        └── temp

    """
    def __init__(self, path=None, new=False):
        if path is not None:
            path = os.path.abspath(path)
        else:
            path = os.getcwd()

        self.config = dict()
        self.path = path
        self.name = os.path.basename(self.path)
        self.index_dir = os.path.join(self.path, index_dir)
        self.config_file = os.path.join(self.path, config_file)
        self.created_file = os.path.join(self.path, created_file)
        self.pid_file = os.path.join(self.path, pid_file)
        self.workflow_file = os.path.join(self.path, workflow_file)
        self.history_dir = os.path.join(self.path, history_dir)
        self.config_dir = os.path.join(self.path, config_dir)
        self.logs_dir = os.path.join(self.path, logs_dir)
        self.results_dir = os.path.join(self.path, results_dir)
        self.temp_dir = os.path.join(self.path, temp_dir)

        if new:
            self.init()

        self.validate()
        self.reload()

        log.info('Loaded project: {}'.format(self.path))

    def __repr__(self):
        return '<Project path={}>'.format(self.path)

    def init(self):
        os.makedirs(self.index_dir, exist_ok=True)
        os.makedirs(self.history_dir, exist_ok=True)
        os.makedirs(self.config_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(self.temp_dir, exist_ok=True)

        if os.path.exists(self.created_file):
            log.info('Reinitialized project: {}'.format(self.path))
        else:
            with open(self.created_file, 'w') as fp:
                info = jetstream.utils.Fingerprint()
                jetstream.utils.yaml_dump(info.serialize(), stream=fp)
            log.info('Initialized project: {}'.format(self.path))

    def validate(self):
        if not os.path.exists(self.path):
            msg = 'Path does not exist: {}'.format(self.path)
            raise NotAProject(msg)

        if not os.path.isdir(self.path):
            raise NotAProject('Not a directory: {}'.format(self.path))

        if not os.path.exists(self.index_dir):
            msg = 'Index dir does not exist {}'.format(self.index_dir)
            raise NotAProject(msg)

        if not os.path.isdir(self.index_dir):
            msg = 'Index dir is not a dir {}'.format(self.index_dir)
            raise NotAProject(msg)

        if not os.path.isdir(self.config_dir):
            msg = 'config/ dir not found at: {}'.format(self.config_dir)
            raise NotAProject(msg)

        if not os.path.isdir(self.history_dir):
            msg = 'history/ dir not found at: {}'.format(self.history_dir)
            raise NotAProject(msg)

        if not os.path.isfile(self.created_file):
            msg = 'Missing: {}'.format(self.created_file)
            raise NotAProject(msg)

    def serialize(self):
        return {k: v for k, v in vars(self).items() if not k.startswith('_')}

    def reload(self):
        """Loads all data files in  <project>/config/ as values in the
        project.config dictionary. """
        if os.path.exists(self.config_file):
            self.config = jetstream.utils.yaml_load(self.config_file)

        if os.path.isdir(self.config_dir):
            config_files = dict()

            for path in jetstream.loadable_files(self.config_dir):
                name = os.path.splitext(os.path.basename(path))[0]
                config_files[name] = jetstream.load_data_file(path)

            self.config.update(config_files)

        with open(self.config_file, 'w') as fp:
            jetstream.utils.yaml_dump(self.config, fp)

    def workflow(self):
        """Load the existing workflow for this project."""
        if os.path.exists(self.workflow_file):
            return jetstream.workflows.load_workflow(self.workflow_file)
        else:
            log.critical('No existing workflow found in this project!')
            return jetstream.Workflow()

    def history(self, paths=False):
        """Find all run ids in this project"""
        runs = []

        for f in os.listdir(self.history_dir):
            if jetstream.run_id_pattern.match(f):
                log.debug('Found run record: {}'.format(f))
                if paths:
                    runs.append(os.path.join(self.history_dir, f))
                else:
                    runs.append(f)
            else:
                log.debug('No run id pattern match, skipping: {}'.format(f))

        return sorted(runs)

    def latest_run(self):
        """Find the latest run folder for this project"""
        try:
            return self.history()[-1]
        except IndexError:
            return None
