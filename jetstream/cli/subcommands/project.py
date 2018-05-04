""" This module contains the cli interface code for the project data utility."""
import sys
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
                        choices=['init', 'config', 'samples'])

    parser.add_argument('args', nargs=argparse.REMAINDER)

    return parser


def init_arg_parser():
    """arg parser for the init action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project init',
    )

    parser.add_argument('path', nargs='?', default='.')

    return parser


def config_arg_parser():
    """arg parser for the data action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project data',
    )

    parser.add_argument('path', nargs='?', default='.')

    parser.add_argument('--format', choices=['yaml', 'json'], default='json')

    parser.add_argument('--json', dest='format',
                        action='store_const', const='json')

    parser.add_argument('--yaml', dest='format',
                        action='store_const', const='yaml')

    parser.add_argument('--pretty', action='store_true', default=False)

    return parser


def init(args=None):
    parser = init_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    jetstream.project.init(args.path)


def config(args=None):
    parser =config_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project(args.path)

    if args.format == 'json':
        if args.pretty:
            print(jetstream.utils.json.dumps(p.config, indent=4))
        else:
            print(jetstream.utils.json.dumps(p.config))
    elif args.format == 'yaml':
        jetstream.utils.yaml.dump(p.config, stream=sys.stdout)


def samples(args=None):
    parser = config_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project(args.path)

    if args.format == 'json':
        if args.pretty:
            print(jetstream.utils.json.dumps(p.samples(), indent=4))
        else:
            print(jetstream.utils.json.dumps(p.samples()))
    elif args.format == 'yaml':
        jetstream.utils.yaml.dump(p.samples(), stream=sys.stdout)


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    if args.action in ('init',):
        init(args=args.args)

    elif args.action in ('config',):
        config(args=args.args)

    elif args.action in ('samples',):
        samples(args=args.args)

    else:
        raise ValueError('Unknown action {}'.format(args.action))
