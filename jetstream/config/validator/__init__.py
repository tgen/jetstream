""" Validators for project config files

These functions are used to validate teh config file key,value pairs against
a schema. Add more value types in the fns module. This should stay relatively
unchanged because the validator fns and schema both load dynamically.
"""

import logging
from collections import Mapping

from . import fns
from .schema import SCHEMA

log = logging.getLogger(__name__)


class SchemaKeyError(Exception):
    """Raised when a key cant be found in schema"""


def get_validator(value_type):
    fn = getattr(fns, 'is_' + value_type.lower())
    return fn


def expected_value_type(key, schema=SCHEMA):
    try:
        res = schema[key]
    except ValueError:
        raise SchemaKeyError(key)
    return res


def check(*args):
    """ check(key, value) or check(mapping)

    Convenience function for validating key-value-pairs against
    the schema """
    if len(args) < 1:
        raise ValueError(check.__doc__)

    elif len(args) == 1:
        mapping = args[0]
        if not isinstance(mapping, Mapping):
            raise ValueError(check.__doc__)
        res = []
        for key, value in mapping.items():
            value_type = expected_value_type(key)
            fn = get_validator(value_type)
            res.append(fn(value))
        return res

    elif len(args) == 2:
        key = args[0]
        value = args[1]
        value_type = expected_value_type(key)
        fn = get_validator(value_type)
        return fn(value)

    else:
        raise ValueError(check.__doc__)

