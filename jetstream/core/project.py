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

# TODO: The idea of storing data about the project in yaml files
# seems like it will be really scalable and powerful, but I am still
# unsure about the best implementation. It is very flexible right now
# but this makes it unclear in some ways. For example, Project.reference()
# will attempt to read a file called "reference.yaml" in the project.
# There is no hard requirement for this file to be present, but some
# tasks in a workflow may expect it. How do we test if all the data
# is present before starting a project? I can only imagine that we need
# to develop strict requirements on the data in a project, but will it
# reduce flexibility?


class MissingProjectData(Exception):
    """Raised when a reference to project data is made for a
    file that does not exist. """
    pass


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

    def __getitem__(self, item):
        """Read-only access to project data files facilitated with
        project['data'] or project.get('data') in a dictionary-like
        pattern. """
        try:
            return self.serialize()[item]
        except KeyError as err:
            raise MissingProjectData(err)

    def get(self, item, fallback=None):
        try:
            self.__getitem__(item)
        except KeyError:
            return fallback

    @property
    def _data_path(self):
        return os.path.join(self.path, 'project.yaml')

    @property
    def _reference_path(self):
        return os.path.join(self.path, 'reference.yaml')

    def _data_files(self):
        """Generator that yields yaml files in the project"""
        for file in os.listdir(self.path):
            if file.endswith(('.yaml', '.yml', '.json', '.tsv', '.csv')):
                yield os.path.join(self.path, file)

    def serialize(self):
        """Returns a dictionary of all data available for this project"""
        data = vars(self)
        for f in self._data_files():
            name = os.path.splitext(os.path.basename(f))[0]
            log.debug('Loading data file: {} as project["{}"]'.format(f, name))

            if f.endswith(('.yaml', '.yml')):
                parsed = utils.yaml_load(f)
            elif f.endswith(('.json',)):
                parsed = utils.json_load(f)
            elif f.endswith(('.csv', '.tsv')):
                parsed = utils.table_to_records(f)
            else:
                log.debug('Skipping unrecognized file type: {}'.format(f))
                continue

            data[name] = parsed

        return data

    def runs(self):
        """Find all run folders for this project"""
        runs = []
        run_data_dir = os.path.join(self.path, RUN_DATA_DIR)
        for i in os.listdir(run_data_dir):
            p = os.path.join(run_data_dir, i)
            if os.path.isdir(p):
                runs.append(i)
        return runs

    def latest_run(self):
        try:
            latest = sorted(self.runs())[-1]
            return latest
        except IndexError:
            return None

    def find(self, item, **kwargs):
        """Returns records in the project[item] if they match on fields given
        as kwargs to this method."""
        return utils.filter_documents(self[item], kwargs)

    def samples(self, **kwargs):
        """ Returns a list of all sample records in project_data. This is
        ephemeral project data generated by joining records from other data
        files, if they are present. It is accomplished by:

        - Starting with any records in under project["samples"]
        - Joining those records with any records found under project["data"]
          key if they match on the "sample_name" property.
        - Filtering the final list of records by any given kwargs

        """
        project_data = self.serialize()
        if 'samples' in project_data:
            samples = {s['sample_name']: s for s in project_data['samples']}
        else:
            samples = {}

        # Sort all data objects from project data into samples
        if 'data' in project_data:
            for data in project_data['data']:
                name = data['sample_name']
                if not name in samples:
                    samples[name] = {
                        'sample_name': name,
                        'data': []
                    }

                samples[name]['data'].append(data)

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
