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

__all__ = [
    "legacy",
    "project",
    "report",
    "pipelines",
    "run",
]
