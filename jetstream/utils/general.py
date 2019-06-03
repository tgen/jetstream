"""Shared utilities"""
import builtins
import confuse
import csv
import fnmatch
import gzip
import importlib
import io
import json
import logging
import os
import subprocess
import sys
import yaml
from collections.abc import Sequence, Mapping
from datetime import datetime
from getpass import getuser
from multiprocessing import cpu_count
from pkg_resources import get_distribution
from socket import gethostname
from uuid import getnode, uuid4

import jetstream

sentinel = object()
log = logging.getLogger(__name__)


class Fingerprint:
    """Generate a snapshot of the system info."""
    def __init__(self, note=None, id=None):
        self.datetime = datetime.utcnow().isoformat()
        self.user = getuser()
        self.version = str(get_distribution("jetstream"))
        self.args = ' '.join(sys.argv)
        self.hostname = gethostname()
        self.pwd = os.getcwd()
        self.note = str(note)
        self.id = id or jetstream.guid()

    def to_dict(self):
        return vars(self)

    def to_yaml(self):
        return yaml_dumps(vars(self))

    def to_json(self):
        return json_dumps(vars(self), sort_keys=True)


class JsonDict(dict):
    """Dict subclass that enforces JSON serializable keys/values.
     This is going to be slower than a dictionary, because it's validating
     every change made to the dict, but will ensure that the items can be
     serialized later."""
    def __init__(self, *args, **kwargs):
        super(JsonDict, self).__init__(*args, **kwargs)
        json_dumps(self)

    def __repr__(self):
        return 'JsonDict({})'.format(super(JsonDict, self).__repr__())

    def __setitem__(self, key, val):
        """When a single item is set, we dont need to revalidate the entire
        object, just the new key/value"""
        json.dumps({key: val})
        return super(JsonDict, self).__setitem__(key, val)

    def update(self, other=None, **kwargs):
        if other is not None:
            json_dumps(other)
            super(JsonDict, self).update(other)
        else:
            json_dumps(kwargs)
            super(JsonDict, self).update(**kwargs)

    def to_dict(self):
        """Return as a standard dictionary"""
        return dict(self)


class Source(str):
    """String subclass that includes a `line_numbers` property for tracking
    the source code line numbers after lines are split up.

    I considered making this an object composed of a string and line number:

    .. code-block:: python

        class Source(object):
            def __init__(self, line_number, data):
              self.line_number = line_number
              self.data = data

        line = Source(0, 'Hello World')


    But, I think it actually complicates most use cases. For example, if the
    source lines were stored in a list, we might want to count a pattern. This
    is easy with a list of strings:

    .. code-block:: python

        res = mylist.count('pattern')


    With a custom class you would be forced to do something like:

    .. code-block:: python

         res = Sum([line for line in lines if line.data == 'pattern'])


    I think this string subclass inheritance pattern is more difficult to
    explain upfront, but it's much is easier to work with downstream. This class
    behaves exactly like a string except in one case: `str.splitlines()` which
    generates a list of Source objects instead of strings. """

    def __new__(cls, data='', line_number=None):
        line = super(Source, cls).__new__(cls, data)
        line.line_number = line_number
        return line

    def splitlines(self, *args, **kwargs):
        """Break source code into a list of source lines."""
        lines = super(Source, self).splitlines(*args, **kwargs)
        lines = [Source(data, line_number=i) for i, data in enumerate(lines)]
        return lines

    def print_ln(self):
        """Format this string along with its line number"""
        return '{}: {}'.format(self.line_number, self)


def coerce_tuple(obj):
    """Coerce an object to a tuple.
    Since strings are sequences, iterating over an object that can be a string
    or a sequence can be difficult. This solves the problem by ensuring scalars
    are converted to sequences. """
    if obj is None:
        return tuple()
    elif isinstance(obj, str):
        return (obj, )
    else:
        return tuple(obj)


def config_stack(*sources):
    """Uses Confuse config flattener to merge together nested config sources.
    The main goal here is to allow multiple sources to merge smoothly rather
    than completely overwrite contents for nested properties.

    In this example, a is the project config and b is the config from command
    arguments. We want the command arguments to overwrite the project config,
    but a standard a.update(b) would not give the desired results

        a = {'foo': {'bar': 42, 'baz': 24}}
        b = {'foo': {'bar': 123}}

        >>> a.update(b)
        >>> a
        {'foo': {'bar': 123}}

        Instead config stack gives us a union:

        >>> config_stack(a, b)
        {'foo': {'bar': 123, 'baz': 24}}


    Warning: this is accomplished with a couple hacks -

    1) To prevent Confuse from searching the filesystem for existing config
    files, an unlikely name is added to the Configuration object when it's
    created.

    2) The resulting config object is json dumped then loaded to get convert
    the ordered dictionaries to standard dicts. So.. high performance is not a
    focus here.

    3) Appending/extending a nested array is still not possible:

        >>> b = {'foo': {'bar': [42, 24,]}}
        >>> a = {'foo': {'bar': [123,]}}
        >>> jetstream.utils.config_stack(a, b)
        {'foo': {'bar': [42, 24]}}

    """
    # Hack to prevent Confuse from accidentally finding a config file on the
    # system when were just trying to create an in-memory config object
    n = 'UNLIKELYAPPNAME-' + str(uuid4())
    conf = confuse.Configuration(n, read=False)
    for s in sources:
        if s:
            if isinstance(s, Mapping):
                conf.set(confuse.ConfigSource(s))
            else:
                err = 'Config sources should return Mapping-type objects'
                raise ValueError(err)

    # Hack to remove the ordereddicts, they're just ugly to look at and when
    # working with primitive data types coming from json/yaml files they're
    # generally useless.
    return json.loads(json.dumps(conf.flatten()))


def dict_update_dot_notation(d, key, value):
    """Updates a keys in dictionaries and allows access to nested dictionaries
    with a dot notation:
    Example:
        >>> d = {'foo': {'bar': 'valA', 'baz': 'valb'}}
        >>> dict_update_dot_notation(d, 'foo.baz', 42)
        >>> d
        {'foo': {'bar': 'valA', 'baz': 42}}
    """
    path = key.split('.')
    k = path.pop(0)

    while 1:
        try:
            nk = path.pop(0)
        except IndexError:
            break

        try:
            d = d[k]
        except KeyError:
            d[k] = dict()
            d = d[k]

        k = nk

    d[k] = value


def dict_lookup_dot_notation(d, path):
    path = path.split('.')

    k = path.pop(0)
    v = d[k]
    d = v

    while path:
        k = path.pop(0)
        v = d[k]
        d = v

    return v


def dynamic_import(path):
    """Imports functions from other libraries when needed"""
    m, _, f = path.rpartition('.')

    try:
        if m and jetstream.settings['dynamic_import'].get():
            mod = importlib.import_module(m)
            return getattr(mod, f)
        else:
            return getattr(builtins, f)
    except AttributeError as e:
        err = f'Failed to import "{path}": Error: {e}'
        raise AttributeError(err) from None


def guess_local_cpus( default=1):
    return cpu_count() or default


def guess_max_forks(default=500):
    """Returns 1/4 of current ulimit -u value. This leaves a good amount of
    room for subprocesses to continue."""
    try:
        res = int(0.25 * int(subprocess.check_output('ulimit -u', shell=True)))
        return res
    except FileNotFoundError as e:
        log.exception(e)
        return default


def is_gzip(path, magic_number=b'\x1f\x8b'):
    """Returns True if the path is gzipped."""
    if os.path.exists(path) and not os.path.isfile(path):
        err = 'This should only be used with regular files because otherwise ' \
              'it will lose some data.'
        raise ValueError(err)

    with open(path, 'rb') as fp:
        if fp.read(2) == magic_number:
            return True
        else:
            return False


def is_multiline(s):
    """ Returns True if the string contains more than one line.
    This tends to perform much better than len(s.splitlines()) or s.count('\n')
    :param s: string to check
    :return: Bool
    """
    start = 0
    lines = 1

    while True:
        idx = s.find('\n', start)

        if idx == -1:
            return lines > 1
        else:
            lines += 1

        if lines > 1:
            return True
        else:
            start = idx + 1


def is_scalar(obj):
    """Returns `True` if the `obj` should be considered a scalar for YAML
    serialization. Strings are an instance of Sequence, so this function might
    look odd, but it does the job."""
    if isinstance(obj, (str,)):
        return True

    if isinstance(obj, (Sequence, Mapping)):
        return False
    else:
        return True


def json_loads(data):
    """Load json data"""
    return json.loads(data)


def json_dumps(obj, *args, sort_keys=True, **kwargs):
    """Attempt to convert `obj` to a JSON string"""
    return json.dumps(obj, sort_keys=sort_keys, *args, **kwargs)


def json_dump(obj, *args, sort_keys=True, **kwargs):
    return json.dump(obj, sort_keys=sort_keys, *args, **kwargs)


def filter_records(records, criteria):
    """Given a list of mapping objects (`records`) and a criteria mapping,
    this function returns a list of objects that match filter criteria."""
    matches = list()
    for i in records:
        for k, v in criteria.items():
            if k not in i:
                log.debug('Drop "{}" due to "{}" not in doc'.format(i, k))
                break
            if i[k] != v:
                log.debug('Drop "{}" due to "{}" not == "{}"'.format(i, k, v))
                break
        else:
            matches.append(i)
    return matches


def find(path, name=None):
    """Similar to BSD find, finds files matching name pattern."""
    for dirname, subdirs, files in os.walk(path):
        for f in files:
            if name is None:
                yield os.path.join(dirname, f)
            elif fnmatch.fnmatch(f, name):
                yield os.path.join(dirname, f)


def load_file(path, filetype=None):
    """Attempts to load a data file from path, raises :ValueError
    if an suitable loader function is not found in get_file_loaders"""
    file_loaders = jetstream.settings['load_file'].get(dict)

    if filetype is not None:
        try:
            loader_name = file_loaders[filetype]
            fn = dynamic_import(loader_name)
            return fn(path)
        except KeyError:
            err = f'No loader defined for {filetype}'
            raise ValueError(err)
    else:
        for ext, loader_name in file_loaders.items():
            if path.endswith(ext):
                fn = dynamic_import(loader_name)
                return fn(path)
        else:
            err = f'No loader extension pattern matched path: {path}'
            raise ValueError(err)


def load_json(path):
    """Load a json file from path"""
    with open(path, 'r') as fp:
        return json.load(fp)


def load_table(path, dialect=None, ordered=False, key=None):
    """Attempts to load a table file in any format. Returns a list of
    dictionaries. This requires the table to have a header row, and does not
    handle comment lines.

    "key" can be used to return a dictionary instead of a list. The key column
    must have a unique value for each row. Here is an
    example of the two different ways a table could be loaded:

    ..
        t = jetstream.utils.load_table('fastqs2.csv')

        [
            {
                'uid': '1',
                'sample': 'sampleA',
                'r1fastq': '/path/to/r1fastq1.gz',
                'r2fastq': '/path/to/r2fastq1.gz'
            },
            {
                'uid': '2',
                'sample': 'sampleA',
                'r1fastq': '/path/to/r1fastq2.gz',
                'r2fastq': '/path/to/r2fastq2.gz'
            },
            {
                'uid': '3',
                'sample': 'sampleA',
                'r1fastq': '/path/to/r1fastq3.gz',
                'r2fastq': '/path/to/r2fastq3.gz'
            }
        ]

        t = jetstream.utils.load_table('fastqs2.csv', key='uid')

        {
            '1': {
                'sample': 'sampleA',
                'r1fastq': '/path/to/r1fastq1.gz',
                'r2fastq': '/path/to/r2fastq1.gz'
                },
            '2': {
                'sample': 'sampleA',
                'r1fastq': '/path/to/r1fastq2.gz',
                'r2fastq': '/path/to/r2fastq2.gz'
                },
            '3': {
                'sample': 'sampleA',
                'r1fastq': '/path/to/r1fastq3.gz',
                'r2fastq': '/path/to/r2fastq3.gz'
                }
        }

    :param path: Path to a table file.
    :param dialect: Instance of csv.Dialect, will be sniffed if dialect is None.
    :param ordered: Return OrderedDict.
    :param key: Return a dictionary instead of a list where items can be
        referenced by their key.
    :returns: If key is given :dict, otherwise a :list
    :rtype: dict, list"""
    with open(path, 'r') as fp:
        data = fp.read()

    if dialect is None:
        dialect = csv.Sniffer().sniff(data)
        log.debug(f'Sniffed dialect: {dialect.__dict__} for "{path}"')

    rows = list(csv.DictReader(data.splitlines(), dialect=dialect))

    if not ordered:
        rows = [dict(r) for r in rows]

    if key:
        new_rows = {}

        for row in rows:
            k = row.pop(key)
            if k in new_rows:
                raise ValueError(f'Duplicate key instance: {k}')
            else:
                new_rows[k] = row

        rows = new_rows

    return rows


def load_yaml(path):
    """Load a yaml file from `path`"""
    with open(path, 'r') as fp:
        return yaml.safe_load(fp.read())


def parse_bool(value):
    """Convert a string value to bool"""
    if value.lower() in ('yes', 'true',):
        return True
    elif value.lower() in ('no', 'false',):
        return False
    else:
        raise TypeError(f'"{value}" could not be interpreted as a boolean')


def read_group(*, ID=None, CN=None, DS=None, DT=None, FO=None, KS=None,
               LB=None, PG=None, PI=None, PL=None, PM=None, PU=None,
               SM=None, strict=True, **unknown):
    """Returns a SAM group header line. This function takes a set of keyword
    arguments that are known tags listed in the SAM specification:

        https://samtools.github.io/hts-specs/

    Unknown tags will raise a `TypeError` unless 'strict' is False

    :param strict: Raise error for unknown read group tags
    :type strict: bool
    :return: str
    """
    if unknown and strict:
        raise TypeError('Unknown read group tags: {}'.format(unknown))

    fields = locals()
    col_order = ('ID', 'CN', 'DS', 'DT', 'FO', 'KS', 'LB', 'PG', 'PI', 'PL',
                 'PM', 'PU', 'SM')

    final = ['@RG']
    for field in col_order:
        value = fields.get(field)
        if value is not None:
            final.append('{}:{}'.format(field, value))

    return '\t'.join(final)


def read_lines_allow_gzip(path):
    """Reads line-separated text files, handles gzipped files and recognizes
    universal newlines. This can cause bytes to be lost when reading from a
    pipe. """
    if is_gzip(path):
        with gzip.open(path, 'rb') as fp:
            data = fp.read().decode('utf-8')
    else:
        with open(path, 'r') as fp:
            data = fp.read()
    lines = data.splitlines()
    return lines


def records_to_csv(records, outpath):
    """Writes records (list of dictionaries) out to a csv file"""
    if os.path.exists(outpath):
        raise FileExistsError(outpath)

    keys = set()
    for record in records:
        keys = keys.union(set(record.keys()))

    log.debug('Found keys: {}'.format(keys))

    with open(outpath, 'w') as fp:
        dw = csv.DictWriter(fp, keys)
        dw.writeheader()
        for record in records:
            # Convert non-scalar values to json strings
            for key, value in record.items():
                if not is_scalar(value):
                    record[key] = json.dumps(value)
        dw.writerows(records)


def remove_prefix(string, prefix):
    if string.startswith(prefix):
        return string[len(prefix):]
    else:
        return string


def yaml_loads(data):
    """Load yaml data"""
    return yaml.safe_load(data)


def yaml_dumps(obj):
    """Attempt to convert `obj` to a YAML string"""
    stream = io.StringIO()
    yaml.dump(obj, stream=stream, default_flow_style=False)
    return stream.getvalue()


def yaml_dump(obj, stream):
    return yaml.dump(obj, stream=stream, default_flow_style=False)
