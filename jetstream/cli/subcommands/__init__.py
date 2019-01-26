"""Subcommands module contains all the argument parsers for Jetstream
commands. When adding to this package, please follow the template
below to give the subcommands some consistent behaviors

Boilerplate for subcommands:

.. code-block:: python

import argparse
import logging

log = logging.getLogger(__name__)


def build_parser():
    parser = argparse.ArgumentParser()
    return parser


def main(args):
    parser = build_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

"""
import sys
import textwrap
import importlib
import pkgutil


def import_submodules(package_name):
    """Import all submodules of a module, recursively

    Thanks to kolypto https://stackoverflow.com/a/25083161/3924113

    :param package_name: Package name
    :type package_name: str
    :rtype: dict[types.ModuleType]
    """
    package = sys.modules[package_name]

    for loader, name, is_pkg in pkgutil.walk_packages(package.__path__):
        yield (name, importlib.import_module(package_name + '.' + name))


cmds = dict(import_submodules(__name__))
__all__ = list(cmds.keys())


def get_module_descriptions(cmds):
    res = list()
    for name, mod in cmds.items():
        desc = {'name': name}
        doc = getattr(mod, '__doc__')
        if doc:
            desc['short_description'] = doc.split('\n\n')[0]
            desc['long_description'] = doc

        res.append(desc)
    return res


def build_summary(width=70):
    title = textwrap.TextWrapper(
        initial_indent=' ' * 2,
        subsequent_indent=' ' * 2,
        width=width
    )

    desc = textwrap.TextWrapper(
        initial_indent=' ' * 4,
        subsequent_indent=' ' * 4,
        width=width
    )

    for mod in get_module_descriptions(cmds):
        name = mod['name']
        description = mod.get('short_description').strip()

        for line in title.wrap(name):
             yield line

        for line in desc.wrap(description or 'No Description'):
            yield line

        yield ''


def summary():
    return '\n'.join([l for l in build_summary()])

