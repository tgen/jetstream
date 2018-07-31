import os
import sys
import io
import csv
import fnmatch
import gzip
import json
import logging
import ulid
import time
import textwrap
import jetstream
from collections.abc import Sequence, Mapping
from datetime import datetime
from getpass import getuser
from socket import gethostname
from uuid import getnode
from urllib.parse import quote as urlquote
from pkg_resources import get_distribution
import yaml


def represent_none(self, _):
    return self.represent_scalar('tag:yaml.org,2002:null', '')


yaml.add_representer(type(None), represent_none)


log = logging.getLogger(__name__)


TEST_RECORDS = [
    {
        'test': 'whitespace',
        'data': ' \t\n\n\x0b\x0c'
    },
    {
        'test': 'punctuation',
        'data': '!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~'
    },
    {
        'test': 'printable',
        'data': '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
                '!"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~ \t\n\n\x0b\x0c'
    },
    {
        'test': 'digits',
        'data': '0123456789'
    },
    {
        'test': 'ascii_letters',
        'data': 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    }
]



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


class JsonDict(dict):
    """ Dict subclass that enforces key/values must be JSON serializable """
    def __init__(self, *args, **kwargs):
        super(JsonDict, self).__init__(*args, **kwargs)

        # If dict was initialized from iterable, recheck the values
        for k, v in self.items():
            self.__setitem__(k, v)

    def __setitem__(self, key, val):
        json.dumps(key)
        json.dumps(val)
        return super(JsonDict, self).__setitem__(key, val)


class LogisticDelay:
    def __init__(self, max=600, inflection=30, sharpness=0.2, ignore=0):
        self.max = max
        self.inflection = inflection
        self.sharpness = sharpness
        self.ignore = ignore
        self.i = 0

    def wait(self):
        self.i += 1
        delay = self.delay(self.i)
        time.sleep(delay)

    def reset(self):
        self.i = 0

    def delay(self, i):
        if i > self.ignore:
            c = self.max
            a = self.inflection
            e = 2.718281828459045
            k = self.sharpness

            return c / (1 + a * e ** (-k * (i - a)))

        else:
            return 0


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


def is_gzip(path, magic_number=b'\x1f\x8b'):
    """Returns True if the path is gzipped."""
    if os.path.exists(path) and not os.path.isfile(path):
        raise OSError("This should only be used with regular files because"
                      "otherwise it will lose some data.")

    with open(path, 'rb') as fp:
        if fp.read(2) == magic_number:
            return True
        else:
            return False


def is_scalar(obj):
    """Returns `True` if the `obj` should be considered
    a scalar for YAML serialization. """
    if isinstance(obj, (str,)):
        return True

    if isinstance(obj, (Sequence, Mapping)):
        return False
    else:
        return True


def remove_prefix(string, prefix):
    if string.startswith(prefix):
        return string[len(prefix):]
    else:
        return string


def json_load(path):
    """Load a json file from path"""
    with open(path, 'r') as fp:
        return json.load(fp)


def json_loads(data):
    """Load json data"""
    return json.loads(data)


def json_dumps(obj, *args, **kwargs):
    """Attempt to convert `obj` to a JSON string"""
    stream = io.StringIO()
    json.dump(obj, fp=stream, sort_keys=True, *args, **kwargs)
    return stream.getvalue()


def json_dump(obj, *args, **kwargs):
    return json.dump(obj, sort_keys=True, *args, **kwargs)


# TODO Handle multi-document yaml files gracefully
def yaml_load(path):
    """Load a yaml file from `path`"""
    with open(path, 'r') as fp:
        return yaml.load(fp)


def yaml_loads(data):
    """Load yaml data"""
    return yaml.load(data)


def yaml_dumps(obj):
    """Attempt to convert `obj` to a YAML string"""
    stream = io.StringIO()
    yaml.dump(obj, stream=stream, default_flow_style=False)
    return stream.getvalue()


def yaml_dump(obj, stream):
    return yaml.dump(obj, stream=stream, default_flow_style=False)


def hash_dict(obj):
    return hash(json_dumps(obj))


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


def table_to_records(path):
    """Attempts to load a table, in any format, as a list of dictionaries

    :param path: Path to a table file
    :return: :list """
    r = list()
    with open(path, 'r') as fp:
        dialect = csv.Sniffer().sniff(fp.readline())
        log.debug('Sniffed delimiter "{}" for "{}"'.format(
            dialect.delimiter, path))
        fp.seek(0)
        reader = csv.DictReader(fp, delimiter=dialect.delimiter)
        for row in reader:
            r.append(dict(row))
    return r


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


def write_test_data(path, dialect='unix'):
    """Writes csv test data out to a file. """
    with open(path, 'w') as fp:
        w = csv.DictWriter(fp, fieldnames=['test', 'data'], dialect=dialect)
        w.writeheader()
        for case in TEST_RECORDS:
            w.writerow(case)
    return path


class Fingerprint(object):
    def __init__(self):
        self.id = str(jetstream.run_id_template.format(ulid.new().str))
        self.datetime = str(datetime.now())
        self.user = str(getuser())
        self.version = str(get_distribution("jetstream"))
        self.sys_version = str(sys.version)
        self.sys_platform = str(sys.platform)
        self.sys_mac = hex(getnode()).upper()
        self.pid = int(os.getpid())
        self.args = ' '.join(sys.argv)
        self.hostname = str(gethostname())
        self.pwd = str(os.getcwd())
        self.parent = str(os.environ.get('JETSTREAM_RUNID'))

    def to_yaml(self):
        return yaml_dumps(vars(self))

    def to_json(self):
        return json.dumps(vars(self))


def fingerprint(to_json=False):
    """Gather system info as a dictionary or JSON string."""
    fp = {
        'id': jetstream.run_id_template.format(ulid.new().str),
        'datetime': str(datetime.now()),
        'user': getuser(),
        'version': str(get_distribution("jetstream")),
        'sys.version': sys.version,
        'sys.platform': sys.platform,
        'sys.mac': hex(getnode()).upper(),
        'pid': os.getpid(),
        'args': sys.argv,
        'hostname': gethostname(),
        'pwd': os.getcwd(),
        'parent': os.environ.get('JETSTREAM_RUNID')
    }

    if to_json:
        return json.dumps(fp)
    else:
        return fp


def find(path, name=None):
    """Similar to BSD find, finds files matching name pattern."""
    for dirname, subdirs, files in os.walk(path):
        for f in files:
            if name is None:
                yield os.path.join(dirname, f)
            elif fnmatch.fnmatch(f, name):
                yield os.path.join(dirname, f)


def cleanse_filename(s):
    return urlquote(s)
