import os
import sys
import stat
import gzip
import logging
import json
from pkg_resources import get_distribution
from socket import gethostname
from getpass import getuser
from uuid import getnode
from datetime import datetime
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
    """ Returns True if the path is gzipped """
    if stat.S_ISFIFO(os.stat(path).st_mode):
        raise OSError("this should not be used with named pipes")

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


# TODO Handle multi-document yaml files gracefully
def yaml_load(path):
    with open(path, 'r') as fp:
        return yaml.safe_load(fp)


def yaml_loads(data):
    return yaml.safe_load(data)


def yaml_dump(obj, path):
    with open(path, 'w') as fp:
        return yaml.dump(obj, stream=fp, default_flow_style=False)


def yaml_dumps(obj):
    return yaml.dump(obj, default_flow_style=False)


def struct(*, action, format, **kwargs):
    dispatcher = {
        'yaml': {
            'load': yaml_load,
            'loads': yaml_loads,
            'dump': yaml_dump,
            'dumps': yaml_dumps
        },
    }

    return dispatcher[format][action](**kwargs)


def fingerprint(to_json=False):
    """Gather system info for recording changes made to projects."""
    fp = {
        'datetime': str(datetime.now()),
        'user': getuser(),
        'version': str(get_distribution(__package__)),
        'sys.version': sys.version,
        'sys.platform': sys.platform,
        'sys.mac': hex(getnode()).upper(),
        'pid': os.getpid(),
        'args': sys.argv,
        'hostname': gethostname(),
        'pwd': os.getcwd()
    }

    if to_json:
        return json.dumps(fp)
    else:
        return fp
