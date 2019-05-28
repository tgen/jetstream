"""Shortcut to "jetstream run" with "--render-only" option"""
import logging
from jetstream.cli.subcommands import run

log = logging.getLogger(__name__)
__doc__ = __doc__+ '\n' + run.__doc__


def arg_parser(parser):
    run.arg_parser(parser)


def main(args):
    log.debug(f'{__name__} {args}')
    args.render_only = True
    return run.main(args)
