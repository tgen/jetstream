"""Command line utility for launching workflows. """
import argparse
import yaml
import json
import logging
from jetstream.core.project import Project
from jetstream.core.workflows.workflow import from_node_link_data
from jetstream.core.run import run_workflow


log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('workflow', help='Path to a workflow file')

    return parser


def launch_single_plugin(plugin, strategy, no_records=False):
    if no_records:
        result = strategy(plugin)
        print(json.dumps(result, indent=4))


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    # Load the project
    p = Project()

    if args.workflow.endswith('.json'):
        with open(args.workflow, 'r') as fp:
            data = json.load(fp)
    else:
        with open(args.workflow, 'r') as fp:
            data = yaml.load(fp.read())

    wf = from_node_link_data(data)

    # Run the workflow
    run_workflow(wf, p)
