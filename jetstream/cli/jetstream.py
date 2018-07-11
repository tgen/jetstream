import argparse
import importlib
import logging
import sys
import traceback
import pkg_resources
import jetstream

log = logging.getLogger()

__version__ = pkg_resources.get_distribution("jetstream").version
log_basic_format = "[\033[4m\033[92m\U0001F335 %(module)10s\033[0m] %(asctime)s: %(message)s"
log_debug_format = "[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s"


def arg_parser():
    main_parser = argparse.ArgumentParser(
        description='Available sub-commands are: {}'.format(get_subcommands()),
        epilog='Use \'jetstream <subcommand> -h/--help\' for help '
               'with specific commands.',
        add_help=False)

    main_parser.add_argument('subcommand', nargs='?', help='subcommand name')

    main_parser.add_argument('-v', '--version', action='version',
                      version=__version__)

    main_parser.add_argument('--verbose', action='store_true')

    main_parser.add_argument('--log-filename', default=None)

    main_parser.add_argument('--log-filemode', default='a')

    main_parser.add_argument('--log-format', default=log_basic_format)

    main_parser.add_argument('--log-level', default='INFO')

    main_parser.add_argument('--ignore-undefined', dest='strict',
                             action='store_false',
                             help='Suppress errors normally raised when a '
                                  'workflow variable is undefined')

    main_parser.add_argument('--template-dir', action='append',
                             help='Specify template search path. If this '
                                  'argument is not set, only built-in templates '
                                  'will be searched. This argument can be used '
                                  'multiple times.')

    main_parser.add_argument('--no-site-templates', dest='site_templates',
                        action='store_false',
                        help='Don\'t include built-in templates in search path')

    return main_parser


def get_subcommands():
    from jetstream.cli.subcommands import __all__ as subcommands
    return ', '.join(subcommands)


def main(args=None):
    parser = arg_parser()
    args, remainder = parser.parse_known_args(args)

    if args.verbose:
        args.log_format = log_debug_format
        args.log_level = 'DEBUG'

    logging.basicConfig(
        filename=args.log_filename,
        filemode=args.log_filemode,
        format=args.log_format,
        level=getattr(logging, args.log_level)
    )

    log.critical('Version {}'.format(__version__))
    log.debug('Cmd args: {}'.format(' '.join(sys.argv)))
    log.debug('{}: {}'.format(__name__, args))

    jetstream.config_environment(
        template_dirs=args.template_dir,
        strict=args.strict,
        include_site_templates=args.site_templates
    )

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
