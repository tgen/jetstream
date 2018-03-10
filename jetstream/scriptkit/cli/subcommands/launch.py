"""Command line utility for launching plugins and workflows"""
import argparse
import json
import logging
import sys
from jetstream.core import Project
from jetstream.workflows import Workflow, launchers

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-p', '--plugin', nargs='+', metavar='PLUGIN_ID',
                        help='Launch plugins from a repo. This takes plugin ids'
                             'in the format: "plugin/path<:revision>"')

    parser.add_argument('-f', '--file', nargs='+',
                        help='Launch plugins directly from files.')

    parser.add_argument('-l', '--launcher',
                        default='default',
                        help=argparse.SUPPRESS,
                        choices=['dry', 'default'])
    return parser


def launch_single_plugin(plugin, strategy, no_records=False):
    if no_records:
        result = strategy(plugin)
        print(json.dumps(result, indent=4))


# TODO Allow arbitrary shell commands to be tracked in job records?
# ex:
# jetstream launch gatk IndelRealigner blah blah
#

def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    launcher = getattr(launchers, args.launcher)
    p = Project()

    if args.file:
        # Generate a new workflow
        wf = Workflow()
        for file_path in args.file:
            plugin_id = 'file://' + file_path
            wf.add_node(plugin_id)

    elif args.plugin:
        # Generate a new workflow
        wf = Workflow()
        for plugin_id in args.plugin:
            wf.add_node(plugin_id)

    else:
        parser.print_help()
        sys.exit(1)

    # Run the workflow
    p.run(wf, launcher=launcher)
