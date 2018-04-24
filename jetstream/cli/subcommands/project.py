""" This module contains the cli interface code for the project data utility."""
import argparse
import logging
import jetstream

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream project',
        description=__doc__
    )

    parser.add_argument('action',
                        choices=['init'])

    parser.add_argument('args', nargs=argparse.REMAINDER)

    return parser


def init_arg_parser():
    """arg parser for the init action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project init',
        description='Add a release'
    )

    parser.add_argument('path', nargs='?', default='.')

    return parser


def init(args=None):
    parser = init_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    jetstream.project.init(args.path)


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    if args.action in ('init',):
        init(args=args.args)

    elif args.action in ('data',):
        # TODO Commandline iterator for samples/data, useful for bash scripts
        raise NotImplementedError

    else:
        raise ValueError('Unknown action {}'.format(args.action))
