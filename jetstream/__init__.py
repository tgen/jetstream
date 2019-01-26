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
from jetstream.projects import Project, NotAProject
from jetstream.runner import Runner
from jetstream.templates import environment, render_template
from jetstream.workflows import Workflow, Task, load_workflow, save_workflow, \
    random_workflow


def context(*, project=None, kvargs=None, separator=None):
    """Returns a dictionary that can be used as the context rendering
    templates with Jinja2"""
    context = {'__project__': project}

    if project:
        config = project.load_config()
        context.update(config)

    if kvargs:
        config = jetstream.kvargs.parse_kvargs(kvargs, separator)
        context.update(config)

    return context


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
    log.debug(f'Logging started: {profile}: {config}')


def represent_str(dumper, data):
    """Allows PyYaml module to dump strings as literals """
    if utils.is_multiline(data):
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
    return dumper.represent_scalar('tag:yaml.org,2002:str', data)

yaml.add_representer(str, representer=represent_str)
