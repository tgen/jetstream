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
    mash='jetstream.cli.subcommands.mash',
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
    use these actions is included in add_config_options_to_parser().

    If this action was assigned to the argument "-c/--config":

    <cmd> -c int:answer 42 -c str:foo 42

    When parsed, these args would add a dictionary to the argparse namespace
    destination "config":

    args = parser.parse_arguments()
    args.config
    {'answer': 42, 'foo': '42'}

    """
    parsers = {}
    _default_parsers = {
        'str': str,
        'int': int,
        'float': float,
        'bool': lambda x: x.lower() == 'true',
        'json': json.loads,
    }

    loaders = {}
    _default_loaders = {}

    def __init__(self, delim=':', *args, **kwargs):
        self.delim = delim
        self.loaders = self.get_loaders()
        self.parsers = self.get_parsers()
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

        if var_type == '':
            var_type = 'str'

        try:
            if var_type.startswith('file:'):
                fn = self.loaders[var_type[5:]]
            else:
                fn = self.parsers[var_type]
        except KeyError:
            msg = f'No function found for "{var_type}". Available types' \
                  f'are: {", ".join(self.get_type_choices())}'
            raise argparse.ArgumentError(self, msg)

        try:
            obj = fn(value)
        except Exception as e:
            tb = traceback.format_exc()
            msg = f'{key} as {var_type}:\n\n{textwrap.indent(tb, "  ")}\n' \
                  f'Error loading config argument: "{key}" with value: ' \
                  f'"{value}", full traceback shown above.'
            raise argparse.ArgumentError(self, msg) from e

        jetstream.utils.dict_update_dot_notation(namespace_dest, key, obj)

    def get_loaders(self):
        loaders = self._default_loaders.copy()
        loaders.update(self.loaders)
        return loaders

    def get_parsers(self):
        parsers = self._default_parsers.copy()
        parsers.update(self.parsers)
        return parsers

    def get_type_choices(self):
        parsers = list(self.parsers.keys())
        loaders = ['file:' + k for k in self.loaders.keys()]
        return parsers + loaders


class UpgradedConfigAction(ConfigAction):
    loaders = jetstream.utils.loaders
    parsers = jetstream.utils.parsers


def add_config_options_to_parser(parser):
    """Adds the -c/--config and -C/--config-file options to an arg parser"""
    act = UpgradedConfigAction(None, [None, None], None)
    config = parser.add_argument_group(
        'config variables',
        description='These options are used to add data that is available for '
                    'rendering templates. Individual items can be added with '
                    '"-c" or data can be loaded in bulk from files with "-C". '
                    'Individual items should follow the syntax '
                    '"-c <[type:]key> <value>". These options can be used '
                    'multiple times. Types can be: '
                    f'{", ".join(act.get_type_choices())}'
    )
    # These options indicate an addition to the args.config object.
    config.add_argument(
        '-c', '--config',
        action=UpgradedConfigAction,
        default={},
        metavar=('TYPE:KEY', 'VALUE'),
        help='add a single template variable'
    )

    # This option allows an monolithic config file to be loaded. Any additional
    # -c config arguments will update this config file.
    config.add_argument(
        '-C', '--config-file',
        metavar='PATH',
        help='load template variables from a config file'
    )

    config.add_argument(
        '--config-file-type',
        help=f'force config file to be loaded as specific file type.'
    )


def arg_parser():
    shared = argparse.ArgumentParser(add_help=False, allow_abbrev=False)

    common = shared.add_argument_group('common options')

    common.add_argument(
        '-l', '--logging',
        choices=jetstream.settings['logging_profiles'].flatten(),
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
        version=jetstream.__version__
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
    log.info(f'Version: {jetstream.__version__}')
    log.debug(f'Command args: {sys.argv}')
    log.debug(f'Settings files:\n{setting_src}')

    # Config item shuffle to get the monolithic config file to serve as the
    # base and then update it with the individual config items. If a config
    # file is loaded that does not yield a dictionary, it will be converted
    # to a dictionary with the contents loaded under __config_file__.
    if hasattr(args, 'config'):
        final = {}
        if args.config_file:
            if args.config_file_type:
                config_file = jetstream.utils.load_file(
                    args.config_file,
                    filetype=args.config_file_type
                )
            else:
                config_file = jetstream.utils.load_file(args.config_file)
            if not isinstance(config_file, dict):
                final = {'__config_file__': config_file}
            else:
                final.update(config_file)
        final.update(args.config)
        args.config = final

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
