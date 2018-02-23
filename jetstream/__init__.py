"""Jetstream is a collection of tools for automating workflows at TGen."""
import pkg_resources
from os import environ

# I'm adding some package wide variables here in order to cleanup the
# namespace of some subpackages.
PLUGIN_DIR = pkg_resources.resource_filename('jetstream', 'plugins/')
PLUGIN_ID_PATTERN = r'(?P<plugin>[^\/]*)\/(?P<path>[^:]*):?(?P<revision>(?<=:)[0-9a-f]{5,40})?$'


# TODO: Move reports.legacy.Project into jetstream.project

# TODO: should project functions walk up the directory tree like git?
# see this https://gist.github.com/zdavkeos/1098474



# TODO Need to discuss most logical organization pattern for these types
# of functions. This is a function that is specific to TGen/Slurm
# infrastructure. Should we break these into a separate module/package
# in order to decouple the engine from infrastructure specifics?
from jetstream.batch_schedulers.slurm import sbatch
from jetstream import config


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

    try:
        project_data = config.load('project.yaml')
        project_name = project_data['name']
    except FileNotFoundError:
        project_name = None

    run_id = environ.get('JETSTREAM_RUN_ID')
    job_name = 'jetstream-{}-{}'.format(project_name, run_id)

    return sbatch(*args, '-J', job_name, stdin_data=final.encode())
