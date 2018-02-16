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


RG_TEMPLATE = 'ID:{ID}\tCN:{CN}\tDS:{DS}\tDT:{DT}\tFO:{FO}\t' \
              'KS:{KS}\tLB:{LB}\tPG:{PG}\tPI:{PI}\tPL:{PL}\t' \
              'PM:{PM}\tPU:{PU}\tSM:{SM}'


def easy_launch(cmd, *args, module_load=None):
    """Launch Slurm jobs with controlled environments via module """
    shebang = '#!/bin/bash\nset -x'

    if module_load:
        final = "{}\nmodule load {} || exit 1\n{}".format(shebang, module_load, cmd)
    else:
        final = "{}\n{}".format(cmd)

    try:
        project_data = config.read('project.yaml')
        project_name = project_data['name']
    except FileNotFoundError:
        project_name = None

    run_id = environ.get('JETSTREAM_RUN_ID')
    job_name = 'jetstream-{}-{}'.format(project_name, run_id)

    return sbatch(*args, '-J', job_name, stdin_data=final.encode())
