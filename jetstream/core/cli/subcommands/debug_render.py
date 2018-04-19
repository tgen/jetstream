"""Render a workflow template with project config data.

*This subcommand is primarily intended for debugging workflow issues. To render,
build, and run a workflow in one step use `jetstream workflow`.*

This must be run inside of a jetstream project. It renders the template with
variables in the project config data, then prints the rendered workflow
template to stdout.

Multiple data files with the same name (basename minus extension) can be used
during the render process in a "last value wins" pattern. Jinja2 will raise e
errors when trying to access nested values for which a parent value is
undefined. However, undefined values in a mapping that exists will only
generate a warning unless "--strict" is used.
"""
import argparse
import logging

import jetstream
from jetstream.core.workflows.builder import render_template

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('template', help='Path to a workflow template')

    parser.add_argument('--strict', action='store_true', default=False)

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project()

    # Load the template
    with open(args.template, 'r') as fp:
        template = fp.read()
    log.critical('Template:\n{}'.format(template))

    # Render the template
    render = render_template(template, obj=p.config, strict=args.strict)
    print(render)
