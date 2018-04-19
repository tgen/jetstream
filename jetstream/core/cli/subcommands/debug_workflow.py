"""Command line utility for running final workflows saved as node-link data.

*This subcommand is primarily intended for debugging workflow issues. To
render, build, and run a workflow in one step use `jetstream workflow`.*
"""
import argparse
import logging

import jetstream
from jetstream.core.project import Project
from jetstream.core.run import run_workflow
from jetstream.core.workflows.workflow import from_node_link_data

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('workflow', help='Path to a workflow file')

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    # Load the project
    p = Project()

    # TODO allow explicit loader declaration
    if args.workflow.endswith('.json'):
        data = jetstream.utils.json_load(args.workflow)
    else:
        data = jetstream.utils.yaml_load(args.workflow)

    log.critical('Workflow data:\n{}'.format(data))
    wf = from_node_link_data(data)

    # Run the workflow
    run_workflow(wf, p)
