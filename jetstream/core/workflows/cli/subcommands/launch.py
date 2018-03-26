"""Command line utility for launching plugins and workflows"""
import argparse
import json
import logging
import sys

from jetstream import utils
from jetstream.core import Project
from jetstream.core.workflows import Workflow, launchers

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    # parser.add_mutually_exclusive_group()
    # This currently causes a bug when using -h here is a discussion
    # on the problem https://bugs.python.org/issue26952

    parser.add_argument('-w', '--workflow', help='Launch a workflow')

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


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    launcher = getattr(launchers, args.launcher)

    # TODO This subcommand could use some work, what if we want to launch
    # plugins from a file AND a repo for example...
    n_opts = sum((
        bool(getattr(args, 'file')),
        bool(getattr(args, 'plugin')),
        bool(getattr(args, 'workflow'))
    ))

    opts_err = 'Select only one option: -f/--file, -p/--plugin, -w/--workflow'

    if n_opts != 1:
        parser.print_help()
        print(opts_err)
        sys.exit(1)

    # Load the project
    p = Project()

    # fh = logging.FileHandler("test.log")
    # fh.setLevel(logging.INFO)
    # log.addHandler(fh)

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

    elif args.workflow:
        # Load the workflow (built beforehand and saved to a file)
        data = utils.yaml_load(args.workflow)
        wf = Workflow.from_node_link_data(data)

    else:
        raise ValueError(opts_err)

    # Run the workflow
    p.run(wf, launcher=launcher)
