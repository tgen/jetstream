import argparse
import traceback
import importlib
import logging
import sys

log = logging.getLogger()

default_log_format = "[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s"


def create_parser():
    main_parser = argparse.ArgumentParser(
        description='Jetstream Toolkit helps streamline scripts. Subcommands '
                    'are: {}'.format(get_subcommands()),
        epilog='Use \'jetstream <subcommand> -h/--help\' for help with a '
               'specific topic.',
        add_help=False)

    main_parser.add_argument('subcommand', nargs='?')

    main_parser.add_argument('-v', '--version', action='store_true')

    main_parser.add_argument('--debug', action='store_true')

    main_parser.add_argument('--log-filename')

    main_parser.add_argument('--log-filemode')

    main_parser.add_argument('--log-format',
                             default=default_log_format)

    main_parser.add_argument('--log-level', default='WARNING')

    return main_parser


def get_subcommands():
    from jetstream.scriptkit.cli.subcommands import __all__ as s
    return s


def main(args=None):
    parser = create_parser()
    args, remaining = parser.parse_known_args(args)

    if args.debug:
        # Alias for lowest logging level
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
        import jetstream
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

        except ModuleNotFoundError as e:
            log.debug(traceback.format_exc())
            parser.print_help()
            log.critical(str(e))
            if args.subcommand != 'help':
                print('Error loading subcommand: {}'.format(args.subcommand))
            sys.exit(1)
