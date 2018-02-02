import os
import shutil
import gzip
import ulid
import logging
from ruamel import yaml

log = logging.getLogger(__name__)

class Source(str):
    """String subclass that includes a "line_numbers" property for tracking
    the source code line numbers after lines are split up.

    I considered making this an object composed of a string and line number:

    ```python
    class Source(object):
        def __init__(self, line_number, data):
          self.line_number = line_number
          self.data = data

    line = Source(0, 'Hello World')
    ```

    But, I think it actually complicates most use cases. For example, if the
    source lines were stored in a list, we might want to count a pattern. This
    is easy with a list of strings:

    ```python
        res = mylist.count('pattern')
    ```

    With a custom class you would be forced to do something like:

    ```python
         res = Sum([line for line in lines if line.data == 'pattern'])
    ```

    I think this string subclass inheritance pattern is more difficult to
    explain upfront, but it's much is easier to work with downstream. This class
    behaves exactly like a string except in one case: str.splitlines() which
    generates a list of Source objects instead of strings. """
    def __new__(cls, data='', line_number=None):
        line = super(Source, cls).__new__(cls, data)
        line.line_number = line_number
        return line

    def splitlines(self, *args, **kwargs):
        lines = super(Source, self).splitlines(*args, **kwargs)
        lines = [Source(data, line_number=i) for i, data in enumerate(lines)]
        return lines

    def print_ln(self):
        return '{}: {}'.format(self.line_number, self)


def read_lines_allow_gzip(path):
    """Reads line-separated text files, handles gzipped files and recognizes
    universal newlines """
    if is_gzip(path):
        with gzip.open(path, 'rb') as fp:
            data = fp.read().decode('utf-8')
    else:
        with open(path, 'r') as fp:
            data = fp.read()
    lines = data.splitlines()
    return lines


def is_gzip(path, magic_number=b'\x1f\x8b'):
    """ Returns True if the path is gzipped """
    with open(path, 'rb') as fp:
        if fp.read(2) == magic_number:
            return True
        else:
            return False


def remove_prefix(string, prefix):
    if string.startswith(prefix):
        return string[len(prefix):]
    else:
        return string


def load_yaml_data(data):
    return yaml.load(data, Loader=yaml.Loader)


def load_yaml(path):
    with open(path, 'r') as fp:
        obj = load_yaml_data(fp.read())
    return obj


def _save_run(run, project):
    target = os.path.join(project, '.jetstream', run['id'])
    lock_file = target + '.lock'
    with open(lock_file, 'w') as fp:
        yaml.dump(project, fp, default_flow_style=False)
    shutil.move(lock_file, target)


def _bind_run_to_project(run, project):
    """ Generates a unique id for the run and saves it to an
    existing project. This will raise a ValueError for runs that
    already have an id, or if project is not a project"""
    if run.get('id') is not None:
        raise ValueError('Run already has id')

    _assert_project(project)

    run['id'] = ulid.new().str
    _save_run(run, project)

    return run['id']


def _is_valid_run_id(id):
    """ Returns True if id is a valid run id """
    try:
        ulid.from_str(id)
        return True
    except (TypeError, ValueError):
        return False


def _is_run(path):
    """ Returns True if path is a valid run """
    if os.path.exists(path) and _is_valid_run_id(path):
        return True
    else:
        return False


def _assert_run(path):
    """ Asserts that path is a run """
    if not _is_run(path):
        raise ValueError('{} is not a jetstream run'.format(path))


def _is_project(path):
    """ Returns True if path is a jetstream project """
    # TODO this should walk up the directory tree like git?
    # see this https://gist.github.com/zdavkeos/1098474
    target = os.path.join(path, '.jetstream/')
    if os.path.exists(target):
        return True
    else:
        return False


def _assert_project(path):
    """ Asserts that path is a jetstream project """
    if not _is_project(path):
        raise ValueError('{} is not a jetstream project'.format(path))


def _find_latest_run(path):
    runs = []
    for d in os.listdir(path):
        if _is_run(d):
            runs.append(d)

    return sorted(runs)[0]


def new_run(path):
    _ = load_project(path)


def load_run(project, run_id):
    """ Loads and returns a specific run from a project """
    _assert_project(project)
    run = os.path.join(project, '.jetstream', run_id)
    _assert_run(run)
    return load_yaml(run)


def load_project(path=None):
    """ Loads and returns the latest run data from a project """
    if path is None:
        path = os.getcwd()

    if not _is_project(path):
        raise ValueError('Not a jetstream project')

    latest = _find_latest_run(path)
    target = os.path.join(path, latest, 'config.yaml')
    data = load_yaml(target)
    return data


def initialize(path=None):
    """ Initialize a Jetstream project in 'path' """
    if path is None:
        path = os.getcwd()
    log.critical('Initializing project in {}'.format(path))
    target = os.path.join(path, '.jetstream/')
    os.makedirs(target, exist_ok=True)
