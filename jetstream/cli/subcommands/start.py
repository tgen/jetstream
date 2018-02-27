"""This module contains the functions that will used to start a new run

A run is:
    One single workflow executed on a project folder


In order to start a run:

    - We need to be inside a project (cwd is a jetstream project)
    - Create a new run. this allows a permanent record of the runs that
      have been performed on a project
    - Load the workflow
    - Start a new run
    - Run requires id, workflow, directory (just for databasing)

To keep things simple for now:

    - there is no path resolution, this application must be started inside
      a valid project directory
    - we do not set up the project directory. that would be performed by the
      application managing the runs.

"""
import json
import argparse
import logging

from jetstream import launchers, workflow, Project

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.set_defaults(action=main)

    parser.add_argument('workflow',
                        help='Path to a workflow file')

    parser.add_argument('-s', '--strategy',
                        default='default')

    return parser

def main(args):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    # Load the strategy (the function that will receive plugins when they're
    # ready to be executed).
    strategy = getattr(launchers, args.strategy)

    # Load the workflow (built beforehand and saved to a file)
    with open(args.workflow, 'r') as fp:
        data = json.load(fp)

    wf = workflow.from_json(data)

    # Start the runner (there may be more than one of these if performance
    # becomes a concern)
    p = Project()
    p.run(wf, strategy=strategy)

