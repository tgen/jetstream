""" This module contains the cli interface code for the config file utility."""
import argparse
import logging
from jetstream import config

log = logging.getLogger(__name__)

def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('action',
                        choices=['convert', 'validate'])

    parser.add_argument('--format', default='json',
                        choices=('json', 'yaml'))

    parser.add_argument('path')

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    if args.action in ('convert',):
        c = config.legacy.load(args.path)
        print(config.serialize(c, format=args.format))

    elif args.action in ('validate',):
        c = config.load(args.path, format=args.format)
        config.validate(c)

    else:
        raise ValueError('Unknown action {}'.format(args.action))