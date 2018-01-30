import os
import sys
import argparse
import jetstream

# TODO add validation to arguments using functions in utils
def validator(value):
    """ Example argument validator function that can be applied with 'type='"""
    if 'condition_is_met':
        return value
    else:
        raise Exception('failed validation')


def add_start_parser(subparser):
    """ Generates the argument parser for the 'start' subcommand. """
    parser = subparser.add_parser('start')
    parser.add_argument('config',
                        help='Path to a config file',
                        type=validator)


def add_update_parser(subparser):
    parser = subparser.add_parser('update')
    parser.add_argument('config',
                        help='Path to a config file',
                        type=validator)


def add_resume_parser(subparser):
    parser = subparser.add_parser('resume')
    parser.add_argument('project',
                        help='Path to a project',
                        type=validator)

def add_config_parser(subparser):
    parser = subparser.add_parser('config')
    parser.add_argument('path',
                        help='Path to a config file')

# TODO add subparsers for plugins module
# Need to be able to sync plugins library, maybe resync,
# or update

# TODO add subparsers for other use cases of the app
# For example, pure shell script plugins can use other
# cli commands to make it easier for accessing project
# data or transforming genomic data

def create_parser():
    main_parser = argparse.ArgumentParser(description=jetstream.__doc__)

    main_parser.add_argument('--always-available')  # Testing precedence for global options

    subparsers = main_parser.add_subparsers(
        title='actions',
        dest='action',
    )

    add_start_parser(subparsers)
    add_update_parser(subparsers)
    add_resume_parser(subparsers)

    return main_parser


def main(args=None):
    parser = create_parser()
    args = parser.parse_args(args)
    if args.action is None:
        parser.print_help()
        print('jetstream -help and jetstream')
    else:
        print(args)


if __name__ == "__main__":
    main()
