import os
import traceback
import logging
import jetstream

log = logging.getLogger(__name__)


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

    Project.list_samples(key=value) returns list of sample records in the that
    has been filtered based on the key/value requirements

    """
    def __init__(self, path=None):
        self.path = path or os.getcwd()
        self.config = dict()
        self.name = os.path.basename(self.path)

        self.workflow_path = os.path.join(self.path, jetstream.project_workflow)
        self.index_path = os.path.join(self.path, jetstream.project_index)
        self.run_history = os.path.join(self.path, jetstream.project_history)
        self.config_path = os.path.join(self.path, jetstream.project_config)
        self.temp_path = os.path.join(self.path, jetstream.project_temp)
        self.log_path = os.path.join(self.path, jetstream.project_logs)
        self.manifest = os.path.join(self.path, jetstream.project_manifest)

        if not os.path.exists(self.path):
            msg = 'Path does not exist: {}'.format(self.path)
            raise jetstream.NotAProject(msg)

        if not os.path.isdir(self.path):
            raise jetstream.NotAProject('Not a directory: {}'.format(self.path))

        if not os.path.exists(self.index_path):
            msg = 'Index dir does not exist {}'.format(self.index_path)
            raise jetstream.NotAProject(msg)

        if not os.path.isdir(self.index_path):
            msg = 'Index dir is not a dir {}'.format(self.index_path)
            raise jetstream.NotAProject(msg)

        self.reload()
        log.info('Loaded project: {}'.format(self.path))

    def __repr__(self):
        return '<Project path={}>'.format(self.path)

    @staticmethod
    def init():
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
            log.info('Initialized project: {}'.format(os.getcwd()))
        else:
            log.info('Reinitialized project: {}'.format(os.getcwd()))

        return Project()

    def serialize(self):
        return {k: v for k, v in vars(self).items() if not k.startswith('_')}

    def reload(self):
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
        for path in jetstream.loadable_files(self.config_path):

            # TODO this might need some more rules in the future
            name = os.path.splitext(os.path.basename(path))[0]

            # Handle legacy configs, see docstring
            if path.endswith('.config') and name == self.name:
                project_legacy_config = path
                continue

            try:
                config[name] = jetstream.load_data_file(path)
            except Exception:
                log.critical('Unable to parse project config: {}'.format(path))
                log.critical(traceback.format_exc())

        if project_legacy_config is not None:
            parsed = jetstream.load_data_file(project_legacy_config)
            config.update(parsed)

        self.config = config

    def pipeline(self):
        """Load the current pipeline for this project. """
        try:
            wf = jetstream.workflows.load_workflow(self.workflow_path)
        except FileNotFoundError:
            wf = jetstream.workflows.Workflow()

        wf.project = self
        return wf

    def runs(self):
        """Find all run folders for this project"""
        runs = []
        run_data_dir = os.path.join(self.path, jetstream.project_history)

        for i in sorted(os.listdir(run_data_dir)):
            try:
                runs.append(Run(self, i))
            except Exception as e:
                log.info('Error loading run: {}'.format(e))

        return runs

    def latest_run(self):
        """Find the latest run folder for this project"""
        try:
            return self.runs()[-1]
        except IndexError:
            return None

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
                try:
                    data_sample_name = record['sample_name']
                except KeyError:
                    msg = 'Unable to link this data record to a sample record: {}'
                    log.warning(msg.format(record))
                    continue

                if data_sample_name not in samples:
                    samples[data_sample_name] = {
                        'sample_name': data_sample_name,
                        'data': []
                    }

                if 'data' not in samples[data_sample_name]:
                    samples[data_sample_name]['data'] = []

                samples[data_sample_name]['data'].append(record)

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

    def add_workflow(self, workflow):
        current = self.pipeline()
        log.info('Current project pipeline {}'.format(current))

        current.compose(workflow)
        log.info('After merging with workflow {}'.format(current))
        return current

    def run(self, workflow, backend=None, *args, **kwargs):
        """ Additional arguments are passed to runner. """
        workflow = self.add_workflow(workflow)
        jetstream.workflows.save(workflow, self.workflow_path)

        if backend is None:
            backend = jetstream.runner.LocalBackend()

        runner = jetstream.AsyncRunner(
            workflow=workflow,
            backend=backend,
            *args, **kwargs)

        rc = runner.start()

        jetstream.workflows.save(workflow, self.workflow_path)

        return rc


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
        path = os.path.join(self.path, filename)

        with open(path, 'w') as fp:
            fp.write(data)
