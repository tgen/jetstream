"""Jetstream CLI

"""
import argparse
import importlib
import json
import logging
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
    draw='jetstream.cli.subcommands.draw',
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
        'json': json.loads
    }

    def __init__(self, delim=':', loaders=None, metavar=None, nargs=2,
                 *args, **kwargs):
        self.delim = delim

        if loaders is None:
            self.loaders = ConfigAction.DEFAULT_LOADERS
        else:
            self.loaders = loaders

        if metavar is None:
            self.metavar = (f'[type{self.delim}]key', 'value')

        if nargs != 2:
            raise ValueError('nargs should always be 2 for ConfigAction')

        super(ConfigAction, self).__init__(*args, nargs=2,  **kwargs)

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
        except Exception:
            tb = traceback.format_exc()
            msg = f'Error loading "{value}" with "{loader}":\n\n' \
                  f'{textwrap.indent(tb, "  ")}\n' \
                  f'Loader exception for "{key}" "{value}", full traceback ' \
                  f'shown above.'
            raise argparse.ArgumentError(self, msg) from None

        jetstream.utils.dict_update_dot_notation(namespace_dest, key, obj)


def add_config_options_to_parser(parser):
    """Adds the -c/--config and -C/--config-file options to an arg parser"""

    # These options indicate an addition to the args.config object.
    parser.add_argument(
        '-c', '--config',
        action=ConfigAction,
        # Paramters seen in ConfigAction.__init__ can also be given here:
        # delim=':'
        help='Add configuration data items individually. These arguments '
             'should follow the syntax "-c <[type:]key> <value>". '
             'This argument can be used multiple times.',
    )

    # This option allows an monolithic config file to be loaded and stored and
    # overwrite the destination args.config. Additionally, the user can supply
    # -c arguments after a -C in order to modify parts of the loaded config
    # file.
    parser.add_argument(
        '-C', '--config-file',
        dest='config',
        type=jetstream.utils.load_file,
        help='load configuration data from a file'
    )


def arg_parser():
    shared = argparse.ArgumentParser(add_help=False, allow_abbrev=False)

    shared.add_argument(
        '-l', '--logging',
        metavar='',
        help='set the logging profile instead of auto-selecting the '
             'logger based on terminal type'
    )

    shared.add_argument(
        '-p', '--project',
        type=load_project
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

    # Dynamically import subcommand modules and add their parsers as
    # subparser to the main parser. This is done just to simplify the
    # process of adding subcommands to the cli. It automatically fills
    # in subparser help text from the docstring of the module.
    subparser = parser.add_subparsers(
        dest='subcommand',
        required=True,
        help=jetstream.cli.subcommands.__doc__
    )

    for cmd, path in _subcommands.items():
        m = importlib.import_module(path)
        p = subparser.add_parser(
            cmd,
            help=m.__doc__.splitlines()[0],
            allow_abbrev=False,
            description=m.__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            parents=[shared,]
        )
        p.set_defaults(func=m.main)
        m.arg_parser(p)

    return parser


def load_project(path=None):
    if path:
        return jetstream.Project(path)
    else:
        try:
            return jetstream.Project()
        except jetstream.ProjectInvalid:
            return None


def main(args=None):
    parser = arg_parser()
    args, remaining = parser.parse_known_args(args)

    jetstream.start_logging(args.logging)
    log.debug(f'Version: {pkg_resources.get_distribution("jetstream")}')
    log.debug(f'sys.argv: {sys.argv}')
    sources = '\n'.join((str(s) for s in jetstream.settings.sources))
    log.debug(f'Settings files:\n{sources}')

    if args.project:
        log.debug(f'Working in {args.project}')
    else:
        log.debug(f'Not working inside of a project!')

    if args.func:
        args.func(args)
