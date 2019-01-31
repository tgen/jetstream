"""Jetstream Command-Line Interface

The jetstream command line tool is actually several sub-commands united
under a single entry point in `jetstream.cli.jetstream.main`. For details on a
sub-command check out the module file in subocommands.
"""
import argparse
import importlib
import jetstream.cli.subcommands
from collections import OrderedDict

# This describes the commands that should be available via this cli. Each
# command should be listed with name as the key and the module import path as
# the value. The OrderedDict will preserve the command order.
_subcommands = OrderedDict(
    build='jetstream.cli.subcommands.build',
    draw='jetstream.cli.subcommands.draw',
    init='jetstream.cli.subcommands.init',
    pipelines='jetstream.cli.subcommands.pipelines',
    project='jetstream.cli.subcommands.project',
    render='jetstream.cli.subcommands.render',
    run='jetstream.cli.subcommands.run',
    settings='jetstream.cli.subcommands.settings',
    tasks='jetstream.cli.subcommands.tasks',
)

[importlib.import_module(p) for c, p in _subcommands.items()]


class Subcommand:
    """Loads a subcommand from subcommands module, handles building the arg
    parser and launching the subcommand.main function. """
    def __init__(self, path, mod, parser):
        self.path = path
        self.mod = mod
        self.parser = parser

    def launch(self, args=None, kvargs=None, project=None):
        args = self.parser.parse_args(args)
        args.kvargs = kvargs
        args.project = project
        self.mod.main(args)


def load_submodule_parsers(parser):
    subparser = parser.add_subparsers(
        dest='subcommand',
        help=jetstream.cli.subcommands.__doc__
    )
    subparser.required = True
    for cmd, path in _subcommands.items():
        m = importlib.import_module(path)
        p = subparser.add_parser(
            cmd,
            help=m.__doc__.splitlines()[0],
            description=m.__doc__,
            formatter_class = argparse.RawDescriptionHelpFormatter,
        )
        p.set_defaults(func=m.main)
        m.arg_parser(p)


def load_submodules():
    modules = {}
    for cmd, path in _subcommands.items():
        mod = importlib.import_module(path)
        p = argparse.ArgumentParser(
            prog=cmd,
            description=mod.__doc__,
            formatter_class = argparse.RawDescriptionHelpFormatter,
        )
        p.set_defaults(func=mod.main)
        mod.arg_parser(p)
        modules[cmd] = Subcommand(path, mod, p)

    return modules