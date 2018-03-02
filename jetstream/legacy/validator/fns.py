"""
Validators search this module for functions to validate
a variable type. Create a new validator with 'is_[type]()'.

Validator functions should return True if the given value
fits the rules for that type, otherwise False.

"""
import os
import re

# Useful regex patterns made available to validators
textField = re.compile('^[a-zA-Z0-9_-]*$')
arrayField = re.compile('^[a-zA-Z0-9_ ,;-]*$')
emailField = re.compile('([^,;\s@]+@[^@]+\.[^,;\s@]+)')
pipelineField = re.compile('(?i)(pegasus|medusa|pecan|chia)')


def is_text(value):
    """ Valid if value contains only alphanumeric, dash, and hyphens """
    return bool(textField.match(value))


def is_array(value):
    """ Valid if value contains only alphanumeric, dash, hyphens, space,
    semicolon, and comma"""
    return bool(arrayField.match(value))


def is_email(value):
    """Valid if value contains at least one @ and one dot, loose email
    address matching. """
    return bool(emailField.match(value))


def is_pipeline(value):
    """Valid True if value is a Pipeline name, case insensitive """
    return bool(pipelineField.match(value))


def is_path(value):
    """Valid if value is a directory that exists"""
    return os.path.isdir(value)


def is_bool(value):
    """ Valid if value is a boolean.
    Lazy booleans are case insensitive matches (partial or whole) to
    truthy (true, yes, 1) or falsey (false, no, 0) values. """
    if isinstance(value, str):
        value = value.lower()
        for b in ('yes', 'true', 'no', 'false', '1', '0'):
            if value in b:
                return True
    elif isinstance(value, int):
        if value in (0, 1):
            return True
    else:
        return False
