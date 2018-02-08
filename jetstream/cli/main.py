import logging
import argparse
from types import ModuleType

from jetstream.cli import subcommands
from jetstream import __doc__ as readme

log = logging.getLogger(__name__)

def create_parser():
    main_parser = argparse.ArgumentParser(description=readme)

    main_parser.add_argument('--log-level', default='WARNING')


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

    logging.basicConfig(
        level=getattr(logging, args.log_level)
    )

    if args.action is None:
        parser.print_help()
        print('jetstream -h to see this message')
    else:
        log.debug(str(args))
        args.action(args)
