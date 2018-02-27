"""Jetstream is a collection of tools for automating workflows at TGen."""
import pkg_resources
import json
from os import environ, getcwd

__version__ = '0.1.0a1'
__author__ = 'Ryan Richholt'
#TODO finish these attributes

PLUGIN_DIR = pkg_resources.resource_filename('jetstream', 'plugins/')
PLUGIN_ID_PATTERN = r'(?P<plugin>[^\/]*)\/(?P<path>[^:]*):?(?P<revision>(?<=:)[0-9a-f]{5,40})?$'

# TODO: should project functions walk up the directory tree like git?
# see this https://gist.github.com/zdavkeos/1098474

from . import workflow, formats, config
from .batch_schedulers import slurm
from .workflow import Workflow
from .project import Project

def read_group(*, ID=None, CN=None, DS=None, DT=None, FO=None, KS=None,
               LB=None, PG=None, PI=None, PL=None, PM=None, PU=None,
               SM=None, strict=True, **unknown):
    """
    Returns a SAM group header line. This function takes a set of keyword
    arguments that are known tags listed in the SAM specification:

        https://samtools.github.io/hts-specs/

    Unknown tags will raise a TypeError unless 'strict' is False

    :param strict: Raise error for unknown read group tags
    :return: Read group string
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


def easy_launch(cmd, *args, module_load=None):
    """Launch Slurm jobs with controlled environments via module """
    shebang = '#!/bin/bash\nset -x'

    if module_load:
        final = "{}\nmodule load {} || exit 1\n{}".format(shebang, module_load, cmd)
    else:
        final = "{}\n{}".format(cmd)

    run_id = environ.get('JETSTREAM_RUN_ID')
    job_name = json.dumps({'jetstream': getcwd(), 'run_id': run_id})

    return slurm.sbatch(*args, '-J', job_name, stdin_data=final.encode())

