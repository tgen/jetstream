import logging
import os
from jetstream import exc, utils
from jetstream.core.settings import profile


log = logging.getLogger(__name__)
RUN_DATA_DIR = profile['RUN_DATA_DIR']

# TODO: should project functions walk up the directory tree like git?
# see this https://gist.github.com/zdavkeos/1098474 but also consider
# this:
# [rrichholt@dback-login1:~]$ git pull
# fatal: Not a git repository (or any parent up to mount point /home)
# Stopping at filesystem boundary (GIT_DISCOVERY_ACROSS_FILESYSTEM not set).


class Project:
    """Internal representation of a project. A project is a directory
    with a run data dir ('project/.jetstream'). Additionally there are
    some data files describing the contents and settings of the project.
    This object provides an interface for easy access to project info.

    Here is a description of the data files that can be present in a
    project their associated getter methods:


    /project.yaml

    Project.meta() will return a dictionary of key:value items in the
    "meta" section of project.yaml

    Project.data() will return a list of all data objects from the "data"
    section of the project.yaml

    Project.data(key=value) will return a list of data objects in the
    "data" section of the project.yaml that has been filtered based on the
    key/value requirements

    Project.samples() will return a list of all samples in the project.yaml
    with their associated data items in an list accessible with sample['data'].
    Note that if sample definitions are included in project.yaml ("samples"
    section is present), they will be joined with the samples_names in data.

    Project.samples(key=value) behaves the same as Project.data() above but
    returns information about samples.


    /reference.yaml

    Project.ref() will return a dictionary of all reference key/values from
    the reference.yaml if it is present in the project directory.

    Project.ref(<key>) will return the value of a single key, identical
    to Project.ref()[<key>].

    """
    def __init__(self, path=None):
        self.path = path or os.getcwd()
        self.name = os.path.basename(self.path)
        self._run_id = ''
        self._run_path = ''

        if not os.path.exists(self.path):
            raise exc.NotAProject('Path does not exist: {}'.format(self.path))

        if not os.path.isdir(self.path):
            raise exc.NotAProject('Not a directory: {}'.format(self.path))

        target = os.path.join(self.path, RUN_DATA_DIR)
        if not os.path.exists(target):
            raise exc.NotAProject('Data dir does not exist {}'.format(target))

        target = os.path.join(self.path, RUN_DATA_DIR)
        if not os.path.isdir(target):
            raise exc.NotAProject('Data dir is not a dir {}'.format(target))

        log.critical('Loaded project {}'.format(self.path))

    def serialize(self):
        return {
            'name': self.name,
            'path': self.path,
            'samples': self.samples(),
            'meta': self.meta(),
            'data': self.data(),
            'ref': self.ref()
        }

    @property
    def _data_path(self):
        return os.path.join(self.path, 'project.yaml')

    @property
    def _reference_path(self):
        return os.path.join(self.path, 'reference.yaml')

    def runs(self):
        run_data_dir = os.path.join(self.path, RUN_DATA_DIR)
        return os.listdir(run_data_dir)

    def latest_run(self):
        try:
            latest = sorted(self.runs())[-1]
            return latest
        except IndexError:
            return None

    def ref(self, *args):
        """Returns a dictionary of all key:values in reference.yaml"""
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
        """Returns a dictionary of all key:values in project.yaml['meta']"""
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
        """Returns a list of all data objects in project.yaml"""
        try:
            project_data = utils.yaml_load(path=self._data_path)
        except FileNotFoundError as err:
            log.warning(err)
            return None

        if kwargs:
            return utils.filter_documents(project_data['data'], kwargs)
        else:
            return project_data['data']

    def samples(self, **kwargs):
        """Returns a list of all sample objects in project.yaml"""
        try:
            project_data = utils.yaml_load(path=self._data_path)
        except FileNotFoundError as err:
            log.warning(err)
            return None

        # Use sample descriptions from project_data if they're available
        if 'samples' in project_data:
            samples = {s['sample_name']: s for s in project_data['samples']}
        else:
            samples = {}

        # Sort all data objects from project data into samples
        for data in self.data():
            sn = data['sample_name']
            if not sn in samples:
                samples[sn] = {
                    'sample_name': sn,
                    'data': []
                }

            samples[sn]['data'].append(data)

        # Turn it into a list so we can filter based on kwargs
        sample_list = list(samples.values())

        if kwargs:
            return utils.filter_documents(sample_list, kwargs)
        else:
            return sample_list


def init():
    if os.path.exists('.jetstream/created'):
        log.critical('Already a project.'.format(os.getcwd()))
    else:
        os.makedirs('.jetstream', exist_ok=True)
        with open('.jetstream/created', 'w') as fp:
            fp.write(utils.yaml.dump(utils.fingerprint()))
        log.critical('Initialized project {}'.format(os.getcwd()))

