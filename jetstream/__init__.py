import logging
import logging.config
import os
import sys

import confuse
import ulid
import yaml


# Configure parallel library dependencies (Used by numpy)
if 'OPENBLAS_NUM_THREADS' not in os.environ:
    os.environ.update(OPENBLAS_NUM_THREADS='1')

if 'MKL_NUM_THREADS' not in os.environ:
    os.environ.update(MKL_NUM_THREADS='1')


# This is an attempt to monkey-patch the confuse.NotFoundError because the
# name of that exception is... confusing. If the class name is updated at
# some point in the future, it would be best to remove this patch.
class ConfigValueNotFound(Exception):
    pass

confuse.NotFoundError = ConfigValueNotFound


# Load settings and configure logging
settings = confuse.LazyConfig('jetstream', __name__)
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())
log.setLevel(1)


# Package module imports
import jetstream.backends
import jetstream.kvargs
import jetstream.pipelines
import jetstream.runner
import jetstream.templates
import jetstream.utils
import jetstream.workflows
from jetstream.projects import Project, ProjectInvalid, new_project
from jetstream.runner import Runner
from jetstream.templates import environment, render_template
from jetstream.workflows import Workflow, Task, load_workflow, save_workflow, \
    random_workflow


def lookup_backend(name=None):
    """Looks up the backend by name or gets the default from the settings. This
    will return the class and also a dictionary of default paramters for
    instantiating the class that can be customized via config file."""
    name = name or jetstream.settings['backend'].get(str)
    params = dict(jetstream.settings['backends'][name].get(dict))
    cls = params.pop('()')
    backend = jetstream.utils.dynamic_import(cls)
    return backend, params


def run_id():
    """Generate a new run ID"""
    return 'js{}'.format(ulid.new().str)


def start_logging(profile=None):
    """Logging is only set up with a NullHandler by default. This function sets
    up logging with the chosen settings profile. """
    if profile is None:
        if sys.stderr.isatty():
            profile = 'interactive'
        else:
            profile = 'default'

    profile = str(profile).lower()
    config = settings['logging_profiles'][profile].get(dict)
    logging.config.dictConfig(config)
    log.debug(f'Logging started: {profile}')


def represent_str(dumper, data):
    """Allows PyYaml module to dump strings as literals """
    if utils.is_multiline(data):
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
    return dumper.represent_scalar('tag:yaml.org,2002:str', data)

yaml.add_representer(str, representer=represent_str)

