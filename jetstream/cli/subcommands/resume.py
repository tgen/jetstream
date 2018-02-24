"""This module contains the functions that will used to resume a run
that crashed without creating a new run object"""
import argparse
import logging

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    return parser


def main(args):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    raise NotImplementedError
