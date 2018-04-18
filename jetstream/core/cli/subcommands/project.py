""" This module contains the cli interface code for the project data utility."""
import argparse
import logging

from jetstream.core import project


log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('action',
                        choices=['data', 'init'])

    parser.add_argument('path', nargs='?')

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    if args.action in ('init',):
        project.init()

    elif args.action in ('data',):
        # TODO Commandline iterator for samples/data, useful for bash scripts
        raise NotImplementedError

    else:
        raise ValueError('Unknown action {}'.format(args.action))
