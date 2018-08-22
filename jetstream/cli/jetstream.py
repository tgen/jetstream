import argparse
import importlib
import logging
import sys
import traceback
import jetstream
from jetstream import log

msg = "%(asctime)s: %(message)s"
color_format = "[\033[4m\033[92m\U0001F335 %(module)10s\033[0m] " + msg
basic_format = "%(module)10s " + msg
verbose_format = "%(name)s.%(funcName)s:%(lineno)d " + msg


def arg_parser():
    main_parser = argparse.ArgumentParser(
        description='Available sub-commands are: {}'.format(get_subcommands()),
        epilog='Use \'jetstream <subcommand> -h/--help\' for help '
               'with specific commands.',
        add_help=False)

    main_parser.add_argument('subcommand', nargs='?', help='subcommand name')

    main_parser.add_argument('-v', '--version', action='version',
                             version=jetstream.__version__)

    main_parser.add_argument('--log-debug', action='store_true',
                             help='Enable debug logging')

    main_parser.add_argument('--log-verbose', action='store_true',
                             help='Enable lowest-level logging')

    main_parser.add_argument('--log-filename', default=None)

    main_parser.add_argument('--log-filemode', default='a')

    main_parser.add_argument('--log-format', default=color_format)

    main_parser.add_argument('--log-level', default='INFO',
                             choices=list(logging._nameToLevel.keys()))

    return main_parser


def get_subcommands():
    from jetstream.cli.subcommands import __all__ as subcommands
    return ', '.join(subcommands)


def main(args=None):
    parser = arg_parser()
    args, remainder = parser.parse_known_args(args)

    log_level = logging._nameToLevel.get(args.log_level, 30)
    log_format = args.log_format
    log_filename = args.log_filename
    log_filemode = args.log_filemode

    if args.log_debug:
        log_level = 10
        log_format = verbose_format

    if args.log_verbose:
        log_level = 1
        log_format = verbose_format

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(log_level)
    stream_handler.setFormatter(logging.Formatter(log_format))
    log.addHandler(stream_handler)

    if log_filename:
        file_handler = logging.FileHandler(log_filename, log_filemode)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(logging.Formatter(basic_format))
        log.addHandler(file_handler)

    log.info('Version {}'.format(jetstream.__version__))
    log.debug('Cmd args: {}'.format(' '.join(sys.argv)))
    log.debug('{}: {}'.format(__name__, args))

    if args.subcommand is None:
        parser.print_help()
    else:
        try:
            # Dynamically import the requested sub-command
            mod = importlib.import_module(
                '.subcommands.' + args.subcommand,
                package=__package__)

            log.debug('Launch {} args: {}'.format(
                mod.__name__, ' '.join(remainder)))

            mod.main(remainder)

        except ModuleNotFoundError:
            log.debug(traceback.format_exc())
            parser.print_help()

            if args.subcommand != 'help':
                print('Error loading subcommand: {}'.format(args.subcommand))

            sys.exit(1)
