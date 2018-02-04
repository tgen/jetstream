import os
import argparse
from types import ModuleType

from jetstream.cli import subcommands
from jetstream import __doc__ as readme


def create_parser():
    main_parser = argparse.ArgumentParser(description=readme)

    main_parser.add_argument('--always-available')  # Testing precedence for global options

    main_parser.add_argument('--cwd', default=os.getcwd(),
                             help='Change the working directory on execution')

    subparsers = main_parser.add_subparsers(
        title='actions',
        dest='action',
    )

    for _, obj in vars(subcommands).items():
        if isinstance(obj, ModuleType):
            obj.arg_parser(subparsers)

    return main_parser


def main(args=None):
    parser = create_parser()
    args = parser.parse_args(args)
    if args.action is None:
        parser.print_help()
        print('jetstream -help and jetstream')
    else:
        args.action(args)
