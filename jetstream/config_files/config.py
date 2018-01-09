import json

from jetstream.components.project import Project
from jetstream import validator

class Source(str):
    """ String subclass that includes a line_numbers property for tracking
    the source code line numbers after lines are split up.

    I considered making this a base object composed of a string and line number:

    ```python
    class SourceLine(object):
        def __init__(self, line_number, data):
          self.line_number = line_number
          self.data = data

    line = SourceLine(0, 'Hello World')
    ```

    But, this actually complicates most use cases. For example, if the source
    lines were stored in a list, we might want to count a pattern. This is easy
    with a list of strings:

    ```python
        res = mylist.count('pattern')
    ```

    With a for a custom class you would be forced to do something like:

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


def read(path, *args, **kwargs):
    """ Read a JSON config file, additional arguments are
    passed to json.loads() """
    with open(path, 'r') as fp:
        data = Source(fp.read())
        parsed = json.loads(data, *args, **kwargs)
        source = data.splitlines()

    project = Project(
        source=source,
        format='json',
        samples=parsed['samples'],
        properties=parsed['properties']
    )

    return project



def validate_metadata(key, value, line=0):
    """ Validate the key-value pair against the schema"""

    is_valid = validator.check(value)

    if not isinstance(is_valid, bool):
        raise Exception('Schema validator error: {}'.format(str(fn)))

    if not is_valid:
        log.warning('Line {}: Schema validation for "{}" failed for "{}"!'.format(
            line, value, fn.__name__))
        # raise SchemaValueError(msg)
    return key, value
