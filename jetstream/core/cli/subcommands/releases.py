"""Command line utility for managing pipeline releases """
import sys
import argparse
import logging
from jetstream.core import releases

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream releases',
        description=__doc__
    )

    parser.add_argument('action',
                        choices=['ls', 'list', 'add', 'remove', 'path'])

    parser.add_argument('args', nargs=argparse.REMAINDER)

    return parser


def add_arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream releases add',
        description='Add a release'
    )

    parser.add_argument('tags', nargs='+', help='Release tags to add')

    return parser



def remove_arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream releases remove',
        description='Remove a release'
    )

    parser.add_argument('tags', nargs='+', help='Release tags to remove')

    return parser


def add(args=None):
    parser = add_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    for tag in args.tags:
        log.critical('Adding {}'.format(tag))
        releases.install_release(tag)


def remove(args=None):
    parser = remove_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    for tag in args.tags:
        log.critical('Removing {}'.format(tag))
        releases.remove(tag)


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    if args.action in ('add',):
        add(args.args)

    elif args.action in ('rm', 'remove'):
        remove(args.args)

    elif args.action in ('path',):
        print(releases.RELEASE_DIR)

    elif args.action in ('ls', 'list'):
        print('\n'.join(releases.list()))

    else:
        parser.print_help()
        log.error('Error! Unrecognized action: {}'.format(args.action))
        sys.exit(1)
