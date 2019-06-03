"""Shortcut to "jetstream run" with "--build-only" option.
"""
import logging
from jetstream.cli.subcommands import run

log = logging.getLogger('jetstream.cli')
__doc__ = __doc__+ '\n' + run.__doc__


def arg_parser(parser):
    run.arg_parser(parser)


def main(args):
    log.debug(f'{__name__} {args}')
    args.build_only = True
    return run.main(args)
