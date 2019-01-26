import sys
import argparse
import pkg_resources
import jetstream
import jetstream.cli.subcommands as subcommands

# Note: This module might make more sense

def arg_parser():
    main_parser = argparse.ArgumentParser(
        description='subcommands:\n\n{}'.format(
            subcommands.summary()),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Use \'jetstream <subcommand> -h/--help\' for help '
               'with specific commands.',
        add_help=False
    )

    main_parser.add_argument(
        'subcommand',
        nargs='?',
        choices=list(subcommands.__all__),
        help=argparse.SUPPRESS
    )

    main_parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=pkg_resources.get_distribution('jetstream').version
    )

    # Global options that configure jetstream.settings. If these are given
    # they take priority over the user/system config file options.
    main_parser.add_argument(
        '--logging',
        metavar='',
        help='set the logging profile instead of auto-selecting the '
             'logger based on terminal type'
    )

    main_parser.add_argument(
        '--kvargs.separator',
        metavar='',
        help='separator for template variable data arguments'
    )

    main_parser.add_argument(
        '--runner.backend',
        metavar='',
        help='set the runner backend used for executing tasks'
    )

    return main_parser


def main(args=None):
    parser = arg_parser()
    args, remainder = parser.parse_known_args(args)

    jetstream.settings.set_args(args, dots=True)
    jetstream.start_logging(args.logging)

    # For backwards compatibility, the --variables argument is treated as a
    # special case and substituted with --file:variables
    try:
        i = remainder.index('--variables')
        remainder[i] = '--file:variables'
    except ValueError:
        pass

    if args.subcommand is None:
        parser.print_help()
        if '-h' in sys.argv or '--help' in sys.argv:
            return
        else:
            raise ValueError('No subcommand given!')
    else:
        mod = getattr(subcommands, args.subcommand)
        mod.main(remainder)
