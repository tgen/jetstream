import sys
import argparse
import logging
from pkg_resources import get_distribution
from jetstream import log, logs, settings
from jetstream.cli import subcommands

description = """Run a Jetstream command.
Available commands are:
{}
""".format(subcommands.summary)


def get_loglevel(value):
    """Determine logging level numeric value from int or level name"""
    try:
        numeric_level = int(value)
    except ValueError:
        numeric_level = logging._nameToLevel[value.upper()]

    return numeric_level


def arg_parser():
    main_parser = argparse.ArgumentParser(
        description='Available commands are:\n\n{}'.format(
            subcommands.summary()),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Use \'jetstream <subcommand> -h/--help\' for help '
               'with specific commands.',
        add_help=False)

    main_parser.add_argument('subcommand', nargs='?',
                             choices=list(subcommands.__all__),
                             help=argparse.SUPPRESS)

    main_parser.add_argument('-v', '--version', action='version',
                             version=get_distribution('jetstream').version)

    main_parser.add_argument('--log-debug', action='store_true',
                             help='Alias for debug log settings')

    main_parser.add_argument('--log-verbose', action='store_true',
                             help='Alias for lowest-level log settings')

    main_parser.add_argument('--log-format', default=None)

    main_parser.add_argument('--log-filename', default=None)

    main_parser.add_argument('--log-filemode', default=None)

    main_parser.add_argument('--log-level', default=None)

    return main_parser


def main(args=None):
    parser = arg_parser()
    args, remainder = parser.parse_known_args(args)

    log_format = args.log_format
    log_level = get_loglevel(args.log_level or settings.get('log_level') or 20)
    log_filename = args.log_filename or settings.get('log_filename')
    log_filemode = args.log_filemode or settings.get('log_filemode')

    if log_filename:
        file_handler = logging.FileHandler(log_filename, log_filemode)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(logging.Formatter(logs.debug_format))
        logs.log.addHandler(file_handler)

    if args.log_debug:
        log_level = 10
        log_format = logs.debug_format

    if args.log_verbose:
        log_level = 1
        log_format = logs.debug_format

    logs.start_logging(format=log_format, level=log_level)
    log.info(f'Version {get_distribution("jetstream")}')
    log.info(f'Cmd args: {" ".join(sys.argv)}')
    log.debug(f'{__name__}: {args}')

    if args.subcommand is None:
        parser.print_help()
        raise ValueError('No subcommand given!')
    else:
        mod = getattr(subcommands, args.subcommand)
        mod.main(remainder)


