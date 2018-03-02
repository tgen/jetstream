from os import environ, getcwd
from jetstream.script_tools.batch_schedulers import slurm


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
    shebang = '#!/bin/bash'

    if module_load:
        final = "{}\nmodule load {} || exit 1\n{}".format(shebang, module_load, cmd)
    else:
        final = "{}\n{}".format(shebang, cmd)

    run_id = environ.get('JETSTREAM_RUN_ID', '')
    job_name = '{}\t{}'.format(getcwd(), run_id)

    return slurm.sbatch(*args, '-J', job_name, stdin_data=final.encode())

