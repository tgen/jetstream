import argparse
import importlib
import logging
import sys
import traceback
import pkg_resources

log = logging.getLogger()

__version__ = pkg_resources.get_distribution("jetstream").version
log_basic_format = "[\033[4m\033[92m\U0001F335 %(module)10s\033[0m] %(asctime)s: %(message)s"
log_debug_format = "[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s"


def create_parser():
    main_parser = argparse.ArgumentParser(
        description='Available subcommands are: {}'.format(get_subcommands()),
        epilog='Use \'jetstream <subcommand> -h/--help\' for help '
               'specific commands.',
        add_help=False)

    main_parser.add_argument('subcommand', nargs='?')

    main_parser.add_argument('-v', '--version', action='version',
                      version=__version__)

    main_parser.add_argument('--debug', action='store_true',
                             help='Alias for lowest level logging')

    main_parser.add_argument('--log-filename')

    main_parser.add_argument('--log-filemode')

    main_parser.add_argument('--log-format', default=log_basic_format)

    main_parser.add_argument('--log-level', default='INFO')

    return main_parser


def get_subcommands():
    from jetstream.cli.subcommands import __all__ as subcommands
    return subcommands


def main(args=None):
    parser = create_parser()
    args, remaining = parser.parse_known_args(args)

    if args.debug:
        args.log_format = log_debug_format
        args.log_level = 'DEBUG'

    logging.basicConfig(
        filename=args.log_filename,
        filemode=args.log_filemode,
        format=args.log_format,
        level=getattr(logging, args.log_level)
    )

    log.debug('Jetstream {}'.format(__version__))
    log.debug('Cmd args: {}'.format(' '.join(sys.argv)))
    log.debug('{}: {}'.format(__name__, args))


    if args.subcommand is None:
        parser.print_help()
        sys.exit(1)

    try:
        # This dynamically imports the requested subcommand
        mod = importlib.import_module(
            '.subcommands.' + args.subcommand,
            package=__package__)

        log.debug('Launch {} args: {}'.format(
            mod.__name__, ' '.join(remaining)))

        mod.main(remaining)

    except ModuleNotFoundError:
        log.debug(traceback.format_exc())
        parser.print_help()
        if args.subcommand != 'help':
            print('Error loading subcommand: {}'.format(args.subcommand))
        sys.exit(1)
