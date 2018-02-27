import os
import argparse
import logging

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.set_defaults(action=main)

    return parser


def main(args):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    proj_dir = os.getcwd()
    os.mkdir(os.path.join(proj_dir, '.jetstream'))

    log.critical('Initialized project {}'.format(proj_dir))
