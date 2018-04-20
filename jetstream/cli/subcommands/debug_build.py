"""Build a workflow from a rendered template

*This subcommand is primarily intended for debugging workflow issues. To render,
build, and run a workflow in one step use `jetstream workflow`.*

If the template has variables they must be rendered before building a workflow.
See `jetstream render`. Templates with Jinja formatting that have not been
rendered will likely throw a yaml parsing exception.
"""
import argparse
import logging

from jetstream import utils
from jetstream.workflows import build_workflow

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('template', help='Path to a workflow template')

    return parser


def main(args=None):
    logging.root.setLevel(logging.DEBUG)
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    # Load the template
    with open(args.template, 'r') as fp:
        template = fp.read()
    log.critical('Template:\n{}'.format(template))

    nodes = utils.yaml_loads(template)
    log.critical('Nodes:\n{}'.format(nodes))

    # Render the template
    render = build_workflow(nodes)
    print(render)
