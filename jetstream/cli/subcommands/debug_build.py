"""Build a workflow from a rendered template file.

Note: This command is primarily intended for debugging template issues. To
render, build, and run a workflow in one step use "jetstream_pipelines".

If the template has variables, they must be rendered before building a workflow.
See "jetstream debug_render".

Templates with Jinja formatting elements that have not been rendered will
(likely) throw a yaml parsing exception. Even if they don't, the result is
probably not what you are intending to create. To be clear, this is low-level
access to the workflow building tools that operate only on rendered workflow
yaml files.
"""
import argparse
import logging

from jetstream import utils
from jetstream.workflows import build_workflow

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('template', help='Path to a workflow template')

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    # Load the template
    with open(args.template, 'r') as fp:
        template = fp.read()
    log.debug('Raw template:\n{}'.format(template))

    nodes = utils.yaml_loads(template)
    log.debug('Node data:\n{}'.format(nodes))

    # Render the template
    render = build_workflow(nodes)
    print(render)
