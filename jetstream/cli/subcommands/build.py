# noinspection PyPackageRequirements
"""Shortcut to "jetstream run" with "--build-only" option.

For complete option listing use "jetstream run -h"

"""
import logging
from jetstream.cli.subcommands import run

log = logging.getLogger('jetstream.cli')


def add_arguments(parser):
    run.add_arguments(parser)


def main(args):
    log.debug(f'{__name__} {args}')
    args.build_only = True
    return run.main(args)
