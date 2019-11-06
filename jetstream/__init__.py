import logging
import logging.config
import os
import sys
import confuse
import ulid
import yaml

__author__ = 'Ryan Richholt'
__email__ = 'rrichholt@tgen.org'
__version__ = '1.666666b1'


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


# Package module imports
from jetstream import backends, pipelines, runner, templates, utils, workflows
from jetstream.projects import Project, init, is_project
from jetstream.runner import Runner
from jetstream.templates import environment, render_template
from jetstream.workflows import Workflow, Task, load_workflow, save_workflow, \
    random_workflow
from jetstream.pipelines import Pipeline, InvalidPipeline, find_pipelines, list_pipelines,\
    get_pipeline


def lookup_backend(name=None):
    """Looks up the backend by name or gets the default from the settings. This
    will return the class and also a dictionary of default paramters for
    instantiating the class that can be customized via config file."""
    name = name or settings['backend'].get(str)
    params = settings['backends'][name].get(dict).copy()
    cls = params.pop('()')
    backend = utils.dynamic_import(cls)
    return backend, params


def guid(formatter=None):
    """Generate a new unique ID"""
    id = ulid.new().str
    if formatter:
        return formatter.format(id=id)
    else:
        return id


def start_logging(profile=None):
    """Logging is only set up with a NullHandler by default. This function sets
    up logging with the chosen settings profile. """
    if profile is None:
        if sys.stdin and sys.stdin.isatty() and sys.stdout and sys.stdout.isatty():
            profile = 'interactive'
        else:
            profile = 'basic'

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

