"""Jetstream CLI

"""
import argparse
import importlib
import json
import logging
import os
import pkg_resources
import sys
import textwrap
import traceback
from collections import OrderedDict
import jetstream
import jetstream.cli.subcommands

log = logging.getLogger('jetstream')

# This describes the commands that should be available via this cli. Each
# command should be listed with name as the key and the module import path as
# the value. The command order will be preserved in the help text.
_subcommands = OrderedDict(
    init='jetstream.cli.subcommands.init',
    run='jetstream.cli.subcommands.run',
    build='jetstream.cli.subcommands.build',
    render='jetstream.cli.subcommands.render',
    tasks='jetstream.cli.subcommands.tasks',
    project='jetstream.cli.subcommands.project',
    pipelines='jetstream.cli.subcommands.pipelines',
    settings='jetstream.cli.subcommands.settings',
)


class ConfigAction(argparse.Action):
    """ConfigAction is an argparse action that allows an arbitrary configuration
    object to be built entirely from command arguments. An example for how to
    use these actions is included in add_config_options_to_parser(). """
    DEFAULT_LOADERS = {
        '': str,
        'str': str,
        'int': int,
        'float': float,
        'bool': lambda x: x.lower() == 'true',
        'json': json.loads,
        'file': jetstream.utils.load_file,
        'file:json': jetstream.utils.load_json,
        'file:yaml': jetstream.utils.load_yaml,
        'file:tsv': jetstream.utils.load_table,
        'file:csv': jetstream.utils.load_table
    }

    def __init__(self, delim=':', loaders=None, nargs=2,
                 *args, **kwargs):
        self.delim = delim

        if loaders is None:
            self.loaders = ConfigAction.DEFAULT_LOADERS
        else:
            self.loaders = loaders

        if nargs != 2:
            raise ValueError('nargs should always be 2 for ConfigAction')

        super(ConfigAction, self).__init__(*args, nargs=2, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None, **kwargs):
        # If this is the first use of the arg, the namespace dest needs
        # to be initialized with an empty dict. Otherwise grab the target
        # out of the namespace
        if getattr(namespace, self.dest, None) is None:
            setattr(namespace, self.dest, dict())
        namespace_dest = getattr(namespace, self.dest)

        key, value = values
        var_type, _, key = key.rpartition(self.delim)

        try:
            loader = self.loaders[var_type]
        except KeyError:
            msg = f'No loader function found for "{var_type}". Available types' \
                  f'are: {",".join(self.loaders.keys())}'
            raise argparse.ArgumentError(self, msg)

        try:
            obj = loader(value)
        except Exception as e:
            tb = traceback.format_exc()
            msg = f'{key} as {var_type}:\n\n{textwrap.indent(tb, "  ")}\n' \
                  f'Loader exception for "{key}" "{value}", full traceback ' \
                  f'shown above.'
            raise argparse.ArgumentError(self, msg) from e

        jetstream.utils.dict_update_dot_notation(namespace_dest, key, obj)

    @staticmethod
    def config_file(value):
        if not os.path.exists(value):
            msg = f'{value} does not exist'
            raise argparse.ArgumentTypeError(msg)

        loader_fn = ConfigAction.DEFAULT_LOADERS['file']

        try:
            return loader_fn(value)
        except Exception as e:
            tb = traceback.format_exc()
            msg = f'\n\n{textwrap.indent(tb, "  ")}\n' \
                  f'Loader exception for "{value}", full traceback ' \
                  f'shown above.'
            raise argparse.ArgumentTypeError(msg) from e


def add_config_options_to_parser(parser):
    """Adds the -c/--config and -C/--config-file options to an arg parser"""
    config = parser.add_argument_group(
        'config variables',
        description='These options are used to add data that is available for '
                    'rendering templates. Individual items can be added with '
                    '"-c" or data can be loaded in bulk from files with "-C". '
                    'Individual items should follow the syntax '
                    '"-c <[type:]key> <value>". These options can be used '
                    'multiple times.'
    )
    # These options indicate an addition to the args.config object.
    config.add_argument(
        '-c', '--config',
        action=ConfigAction,
        metavar=('TYPE:KEY', 'VALUE'),
        help='add a single template variable'
    )

    # This option allows an monolithic config file to be loaded and stored and
    # overwrite the destination args.config. Additionally, the user can supply
    # -c arguments after a -C in order to modify parts of the loaded config
    # file.
    config.add_argument(
        '-C', '--config-file',
        dest='config',
        metavar='PATH',
        type=ConfigAction.config_file,
        help='load template variables from a file'
    )


def arg_parser():
    shared = argparse.ArgumentParser(add_help=False, allow_abbrev=False)

    common = shared.add_argument_group('common options')

    common.add_argument(
        '-l', '--logging',
        choices=jetstream.settings['logging_profiles'].get(dict),
        help='set the logging profile [interactive]'
    )

    common.add_argument(
        '-p', '--project',
        help='use project directory [set automatically if current directory is '
             'a jetstream project]'
    )

    common.add_argument(
        '--pipeline',
        help='use pipeline directory [set automatically if using pipelines command]'
    )

    parser = argparse.ArgumentParser(
        prog='jetstream',
        allow_abbrev=False,
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[shared,],
        epilog='Use \'jetstream <subcommand> -h/--help\' for help '
               'with specific commands.',
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version=pkg_resources.get_distribution('jetstream').version
    )

    parser.set_defaults(func=None)

    # Dynamically import subcommand modules and add subparser to the
    # main parser for each subcommand. This is done just to simplify the
    # process of adding new subcommands to the cli. It automatically fills
    # in subparser help text from the docstring of the module, and calls
    # the main function with the parsed Namespace.
    subparser = parser.add_subparsers(
        dest='subcommand',
        help=jetstream.cli.subcommands.__doc__
    )

    for cmd, path in _subcommands.items():
        m = importlib.import_module(path)
        p = subparser.add_parser(
            cmd,
            help=m.__doc__.splitlines()[0],
            allow_abbrev=False,
            description=m.__doc__,
            parents=[shared,]
        )
        p.set_defaults(func=m.main)
        m.add_arguments(p)

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)

    jetstream.start_logging(args.logging)
    setting_src = '\n'.join((str(s) for s in jetstream.settings.sources))
    log.info(f'Version: {pkg_resources.get_distribution("jetstream")}')
    log.debug(f'Command args: {sys.argv}')
    log.debug(f'Settings files:\n{setting_src}')

    if args.project:
        args.project = jetstream.Project(args.project)
    else:
        cwd = os.getcwd()
        if jetstream.projects.is_project(cwd):
            args.project = jetstream.Project(cwd)

    if args.pipeline:
        args.pipeline = jetstream.Pipeline(args.pipeline)

    if args.func:
        args.func(args)
    else:
        parser.error('the following arguments are required: subcommand')
