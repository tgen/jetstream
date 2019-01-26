"""Create a new project or reinitialize an existing project.

This command can be used to create a new Jetstrema project with the
recommended project folders. If config files are included, they will
be copied into the project config directory."""
import os
import argparse
import jetstream
from jetstream import log


def arg_parser():
    """Argument parser for the init action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project init',
        description=__doc__
    )

    parser.add_argument(
        'path',
        nargs='?',
        default=os.getcwd(),
        help='Path to a Jetstream project'
    )

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    os.makedirs(args.path, exist_ok=True)
    p = jetstream.Project(path=args.path)
    p.create()

    log.info('Done')
