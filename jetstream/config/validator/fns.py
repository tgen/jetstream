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
    """ True if value contains only alphanumeric, dash, and hyphens """
    return bool(textField.match(value))


def is_array(value):
    """ Returns true if value contains only alphanumeric, dash, hyphens, space, semicolon, and comma"""
    return bool(arrayField.match(value))


def is_email(value):
    """ True if value contains at one @ and one dot, loose email address matching"""
    return bool(emailField.match(value))


def is_pipeline(value):
    """ Returns True if value is a Pipeline name, case insensitive """
    return bool(pipelineField.match(value))


def is_path(value):
    """ Returns True if value is a directory """
    return os.path.isdir(value)


def is_bool(value):
    """ Returns True if value is a valid lazy boolean.
    Lazy booleans are case insensitive matches (partial or whole) to
    truthy or falsey values: true, false, yes, no, 1 or 0. """
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
