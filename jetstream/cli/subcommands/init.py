import os
import shutil
import argparse
import jetstream
from jetstream import log


def init_arg_parser():
    """Argument parser for the init action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project init',
        description='Initialize a Jetstream project'
    )

    parser.add_argument('path', nargs='?', default=os.getcwd(),
                        help='Path to a Jetstream project')

    parser.add_argument('-c', '--config', nargs='*', default=list(),
                        help='Copy files into the new \'<project>/config/\'')

    return parser


def main(args=None):
    parser = init_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    os.makedirs(args.path, exist_ok=True)
    p = jetstream.Project(path=args.path, new=True)

    for file in args.config:
        log.info('Copying: {} -> {}'.format(file, p.config_dir))
        shutil.copy(file, p.config_dir)
