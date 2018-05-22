import os
import traceback
import logging
import jetstream
from jetstream.runners import BaseRunner

log = logging.getLogger(__name__)


class NotAProject(Exception):
    pass


class NotARun(Exception):
    pass


class Project:
    """Interact with a Jetstream project.

    Internal representation of a project. A project is a directory
    with a run data dir ('project/.jetstream'). Additionally, there are
    some data files describing the contents and settings of the project.
    This object provides an interface for easy access to project info.

    Here is a description of the data files that can be present in a
    project their associated getter methods:

    Project.samples() will return a list of all sample records in the project
    data with their associated data records in an list accessible with
    sample['data']. Note that if sample definitions are included in
    project.data['samples'], they will be joined with the samples_names in
    data records.

    Project.samples(key=value) returns list of sample records in the that has
    been filtered based on the key/value requirements

    """
    def __init__(self, path=None):
        self.path = path or os.getcwd()
        self.config = dict()
        self.name = os.path.basename(self.path)

        self.pid_file = os.path.join(self.path, jetstream.project_pid_file)
        self.workflow_path = os.path.join(self.path, jetstream.project_workflow)
        self.index_path = os.path.join(self.path, jetstream.project_index)
        self.run_history = os.path.join(self.path, jetstream.project_history)
        self.config_path = os.path.join(self.path, jetstream.project_config)
        self.temp_path = os.path.join(self.path, jetstream.project_temp)
        self.log_path = os.path.join(self.path, jetstream.project_logs)
        self.manifest = os.path.join(self.path, jetstream.project_manifest)

        if not os.path.exists(self.path):
            raise NotAProject('Path does not exist: {}'.format(self.path))

        if not os.path.isdir(self.path):
            raise NotAProject('Not a directory: {}'.format(self.path))

        if not os.path.exists(self.index_path):
            raise NotAProject('Index dir does not exist {}'.format(
                self.index_path))

        if not os.path.isdir(self.index_path):
            raise NotAProject('Index dir is not a dir {}'.format(
                self.index_path))

        self._load_project_config_files()
        log.critical('Loaded project: {}'.format(self.path))

    def __repr__(self):
        return '<Project path={}>'.format(self.path)

    @staticmethod
    def init(path=None):
        cwd = os.getcwd()
        try:
            if path is not None:
                os.makedirs(path, exist_ok=True)
                os.chdir(path)

            os.makedirs(jetstream.project_index, exist_ok=True)
            os.makedirs(jetstream.project_history, exist_ok=True)
            os.makedirs(jetstream.project_config, exist_ok=True)
            os.makedirs(jetstream.project_temp, exist_ok=True)
            os.makedirs(jetstream.project_logs, exist_ok=True)

            created_path = os.path.join(jetstream.project_index, 'created')
            if not os.path.exists(created_path):
                with open(created_path, 'w') as fp:
                    created = jetstream.utils.fingerprint()
                    jetstream.utils.yaml.dump(created, stream=fp)
                log.critical('Initialized project: {}'.format(os.getcwd()))
            else:
                log.critical('Reinitialized project: {}'.format(os.getcwd()))

        finally:
            os.chdir(cwd)

    def serialize(self):
        return {k: v for k, v in vars(self).items() if not k.startswith('_')}

    def _load_project_config_files(self):
        """Loads all data files in the project/config as values in the
        project.config dictionary.

        Legacy configs are handled differently than other data files. If the
        name of the config file matches the name of the project, the values
        for "meta" and "data" are added directly to project.data and will
        overwrite any previous values. Legacy config files with names other
        than the current project name are handled just like other data files.
        """
        # TODO Smarter merging of config files when some records can
        # be stored in multiple files https://pypi.org/project/jsonmerge/
        config = dict()
        project_legacy_config = None
        for path in loadable_files(self.config_path):
            name = path_root(path)

            # Handle legacy configs, see docstring
            if path.endswith('.config') and name == self.name:
                project_legacy_config = path
                continue

            try:
                config[name] = load_data_file(path)
            except Exception:
                log.warning('Unable to parse project config: {}'.format(path))
                log.debug(traceback.format_exc())

        if project_legacy_config is not None:
            parsed = load_data_file(project_legacy_config)
            config.update(parsed)

        self.config = config

    def load_workflow(self):
        try:
            wf = jetstream.workflows.load_workflow(self.workflow_path)
        except FileNotFoundError:
            wf = jetstream.workflows.Workflow()

        return wf

    def runs(self):
        """Find all run folders for this project"""
        runs = []
        run_data_dir = os.path.join(self.path, jetstream.project_history)

        for i in sorted(os.listdir(run_data_dir)):
            try:
                runs.append(Run(self, i))
            except Exception as e:
                log.critical('Error loading run: {}'.format(e))

        return runs

    def latest_run(self):
        try:
            return self.runs()[-1]
        except IndexError:
            return None

    def new_run(self):
        """Create a new Run"""
        data = jetstream.utils.fingerprint()
        run_id = data['id']
        path = os.path.join(self.path, jetstream.project_history, run_id)

        os.makedirs(path, exist_ok=True)

        with open(os.path.join(path, 'created'), 'w') as fp:
            jetstream.utils.yaml.dump(data, fp)

        return Run(self, run_id)

    def samples(self):
        """Build a dictionary of all sample records in project.config.
       This is ephemeral project data generated by joining records from several
       project config files (if they are present):

       - Starting with any records in under project.config["samples"]
       - Joining those records with any found under project.config["data"]
         key if they match on the "sample_name" property.
       - Filtering the results by any given kwargs

       """
        if 'samples' in self.config:
            samples = {s['sample_name']: s for s in self.config['samples']}
        else:
            samples = {}

        # Sort all data objects from project.config['data'] into samples
        if 'data' in self.config:
            for record in self.config['data']:
                record_sample_name = record['sample_name']

                if record_sample_name not in samples:
                    samples[record_sample_name] = {
                        'sample_name': record_sample_name,
                        'data': []
                    }

                if 'data' not in samples[record_sample_name]:
                    samples[record_sample_name]['data'] = []

                samples[record_sample_name]['data'].append(record)

        return samples

    def list_samples(self, **kwargs):
        """Returns a (filtered) list of all sample records in project.config.

        This list contains the values of all records found by Project.samples.
        It will be filtered to contain only records where attributes match the
        given kwargs. """
        sample_list = list(self.samples().values())

        if kwargs:
            return jetstream.utils.filter_records(sample_list, kwargs)
        else:
            return sample_list

    def render(self, template, additional_data=None):
        if additional_data is None:
            additional_data = dict()

        temp = jetstream.env.get_template_with_source(template)
        return temp.render(project=self, **additional_data)

    def run(self, template, additional_data=None, runner_class=BaseRunner):
        if additional_data is None:
            additional_data = dict()

        # The changes to workflow composition caused nested workflows
        # to generate endless loops. Disabling nested workflows via
        # pid files for now until I can rework the workflow composition.
        if os.path.exists(self.pid_file):
            raise FileExistsError(self.pid_file)
        else:
            with open(self.pid_file, 'w') as fp:
                jetstream.utils.yaml.dump(jetstream.utils.fingerprint(), fp)

        try:
            run = self.new_run()
            temp = jetstream.env.get_template_with_source(template)
            run.save(jetstream.utils.yaml_dumps(temp.source), 'template')

            tasks = temp.render(project=self, **additional_data)
            run.save(tasks, 'tasks')

            workflow = jetstream.workflows.build_workflow(tasks)
            run.save(str(workflow), 'workflow')

            existing_workflow = self.load_workflow()
            existing_workflow.compose(workflow)
            existing_workflow.project = self
            existing_workflow.auto_save = True

            existing_workflow.retry()

            runner = runner_class(existing_workflow)
            return runner.start()

        finally:
            os.remove(self.pid_file)


class Run(object):
    def __init__(self, project, run_id):
        self.id = run_id
        self.project = project
        self.path = os.path.join(project.path, jetstream.project_history, self.id)
        self.info = self._load_command_info()

    def __repr__(self):
        return '<Run {}: {}>'.format(self.id, self.info.get('datetime'))

    def _load_command_info(self):
        return jetstream.utils.yaml_load(os.path.join(self.path, 'created'))

    def save(self, data, filename):
        """Save something to the run history"""
        with open(os.path.join(self.path, filename), 'w') as fp:
            fp.write(data)


def init(path=None):
    return Project.init(path)


def loadable_files(directory):
    """Generator yields all files we can load (see data_loaders) """
    for file in os.listdir(directory):
        path = os.path.join(directory, file)
        if os.path.isfile(path) \
                and path.endswith(tuple(jetstream.data_loaders.keys())):
            yield path


def path_root(path):
    """Returns path minus directories and extension"""
    # TODO this might need some more rules in the future
    return os.path.splitext(os.path.basename(path))[0]


def load_data_file(path):
    """Attempts to load a data file from path, raises :ValueError
    if an suitable loader function is not found in data_loaders"""
    for ext, fn in jetstream.data_loaders.items():
        if path.endswith(ext):
            loader = fn
            break
    else:
        raise ValueError('No loader fn found for {}'.format(path))

    log.debug('Loading {} with {}.{}'.format(
        path, loader.__module__, loader.__name__))

    return loader(path)
