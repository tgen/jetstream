"""Command line utility for one-step render and run of workflow templates"""
import argparse
import logging

from jetstream import utils
from jetstream.core.project import Project
from jetstream.core.run import run_workflow
from jetstream.core.workflows.builder import render_template, build_workflow

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

    # Load the project, ensure we're working in a valid project
    p = Project()

    # Load the template
    log.critical('Loading workflow template... {}'.format(args.template))
    with open(args.template, 'r') as fp:
        template = fp.read()
    log.debug('Template:\n{}'.format(template))

    # Render variables in the workflow template
    log.critical('Rendering template... strict: {}'.format(args.strict))
    rendered_template = render_template(
        template=template,
        obj=p.config,
        strict=args.strict
    )

    # Rendered template is a yaml format array of nodes
    # we load this in with the yaml library, then build a
    # workflow from the nodes
    log.critical('Loading node-link data...')
    nodes = utils.yaml_loads(rendered_template)

    log.critical('Building workflow...')
    wf = build_workflow(nodes)

    log.critical('Starting run...')
    # Now we run the workflow in the project
    run_workflow(wf, p)
