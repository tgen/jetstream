import sys
import logging
import argparse
import importlib

import jetstream

log = logging.getLogger(__name__)

SUBCOMMANDS = (
    'config',
    'launch',
    'plugins',
    'report',
    'resume',
    'start',
    'help'
)

def create_parser():
    main_parser = argparse.ArgumentParser(
        description=jetstream.__doc__,
        epilog='Use \'jetstream <subcommand> -h/--help\' for help with a '
               'specific topic.',
        add_help=False)

    main_parser.add_argument('subcommand', nargs='?', choices=SUBCOMMANDS)

    main_parser.add_argument('-v', '--version', action='store_true')

    main_parser.add_argument('--debug', action='store_true')

    main_parser.add_argument('--log-filename')

    main_parser.add_argument('--log-filemode')

    main_parser.add_argument('--log-format')

    main_parser.add_argument('--log-level', default='WARNING')

    return main_parser


def main(args=None):
    parser = create_parser()
    args, remaining = parser.parse_known_args(args)

    if args.debug:
        # Alias for lowest level logging
        args.log_level = 'DEBUG'

    logging.basicConfig(
        filename=args.log_filename,
        filemode=args.log_filemode,
        format=args.log_format,
        level=getattr(logging, args.log_level)
    )

    log.debug(sys.argv)
    log.debug('{}: {}'.format(__name__, args))

    if args.version:
        print(jetstream.__version__)
        sys.exit(0)
    else:
        if args.subcommand is None:
            parser.print_help()
            sys.exit(1)

        try:
            # This dynamically imports the requested subcommand
            mod = importlib.import_module(
                '.subcommands.' + args.subcommand,
                package=__package__)

            log.debug('Launch {} remaining args: {}'.format(mod, remaining))
            mod.main(remaining)

        except ModuleNotFoundError:
            parser.print_help()
            if args.subcommand != 'help':
                print('Error! Unknown subcommand: {}'.format(args.subcommand))
            sys.exit(1)
