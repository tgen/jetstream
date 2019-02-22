"""Jetstream pipeline development toolkit

TODO: all arguments following '--' will be ignored by the argument parser
and used as template rendering data.
"""
import argparse
import logging
import os
import pkg_resources
import sys
import jetstream as js  # The jetstream.py file here conflicts with package name
import jetstream.cli as cli

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream',
        allow_abbrev=False,
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Use \'jetstream <subcommand> -h/--help\' for help '
               'with specific commands.',
        add_help=False,
    )

    parser.add_argument(
        'subcommand',
        nargs='?',
        choices=list(cli._subcommands.keys())
    )

    parser.add_argument(
        '-k',
        '--kvargs.separator',
        metavar='',
        help='separator for template variable data arguments'
    )

    parser.add_argument(
        '-l',
        '--logging',
        metavar='',
        help='set the logging profile instead of auto-selecting the '
             'logger based on terminal type'
    )

    parser.add_argument(
        '-p',
        '--project',
        help='If the cwd is a project, it will be loaded automatically. '
             'Otherwise, a path to a project can be specified, and the run '
             'will start in that project directory.'
    )

    parser.add_argument(
        '-r',
        '--runner.backend',
        metavar='',
        help='set the runner backend used for executing tasks'
    )

    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=pkg_resources.get_distribution('jetstream').version
    )

    parser.add_argument(
        '-h',
        '--help',
        action='store_true'
    )

    return parser


def main(args=None):
    # This forces any arguments following "--" to be completely ignored by
    # the argparsers so that we can use them for template context
    args = args or sys.argv[1:]
    try:
        i = args.index('--')
        kvargs = args[i + 1:]
        args = args[:i]
    except ValueError:
        kvargs = []

    parser = arg_parser()
    args, remaining = parser.parse_known_args(args)
    kvargs = js.kvargs.parse_kvargs(kvargs)

    if args.help:
        if args.subcommand:
            remaining.insert(0, '-h')
        else:
            return parser.print_help()

    if args.subcommand is None:
        parser.error('the following arguments are required: subcommand')

    js.settings.set_args(args, dots=True)
    js.start_logging(args.logging)

    log.info(f'{pkg_resources.get_distribution("jetstream")}')
    log.debug(f'sys.argv: {sys.argv}')
    log.debug(f'config sources: {js.settings.sources}')

    if args.project:
        project = js.Project(args.project)
    else:
        try:
            project = js.Project()
        except js.ProjectInvalid:
            project = None

    if project:
        log.info(f'Working in {project}')
    else:
        log.info(f'Not working inside of a project!')

    mod = cli.load_submodules()[args.subcommand]
    mod.launch(remaining, kvargs, project)
