""" Validators for project config files

These functions are used to validate teh config file key,value pairs against
a schema. Add more value types in the fns module. This should stay relatively
unchanged because the validator fns and schema both load dynamically.
"""

import logging
from collections import Mapping

from . import fns
from .schema import SCHEMA

# SCHEMA keys are case insensitive
SCHEMA = {k.lower():v for k,v in SCHEMA.items()}

log = logging.getLogger(__name__)


class KeyNotInSchema(Exception):
    """Raised when a key cant be found in schema"""

class FailedValidation(Exception):
    """Raised when a key-value pair failed validation"""


def get_validator(value_type):
    fn = getattr(fns, 'is_' + value_type.lower())
    return fn


def expected_value_type(key, schema=SCHEMA):
    try:
        return schema[key.lower()]
    except KeyError:
        raise KeyNotInSchema(key) from None


def check_kv(key, value):
    """
    :param key: Key to check
    :param value: Value to validate
    :return: True if validator function returns true
    """
    value_type = expected_value_type(key)
    fn = get_validator(value_type)

    if not fn(value):
        msg = 'Failed Validator\n' \
              'Key: {}\n' \
              'Value: {}\n' \
              'Type: {}\n' \
              'Vaidator: {}\n' \
              'Help: {}'.format(
            key, value, value_type, str(fn), fn.__doc__)
        raise FailedValidation(msg)
    else:
        return True


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

        good = True
        for key, value in mapping.items():
            try:
                check_kv(key, value)
            except (KeyNotInSchema, FailedValidation) as e:
                log.warning(e)
                good = False

        if not good:
            raise FailedValidation


    elif len(args) == 2:
        check_kv(*args)

    else:
        raise ValueError(check.__doc__)

