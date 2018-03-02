""" This module contains the cli interface code for the project data utility."""
import argparse
import logging
from jetstream import utils
import jetstream.projects as projects
from jetstream.legacy import config

log = logging.getLogger(__name__)

def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('action',
                        choices=['data', 'init', 'legacy'])

    parser.add_argument('path', nargs='?')

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    if args.action in ('init',):
        projects.init()

    elif args.action in ('legacy',):
        c = config.load(args.path)
        print(utils.yaml_dumps(c))

    elif args.action in ('data',):
        raise NotImplementedError

    else:
        raise ValueError('Unknown action {}'.format(args.action))