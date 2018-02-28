"""Command line utility for launching plugins and workflows"""
import sys
import argparse
import json
import logging
from jetstream import launchers, workflow, Project, utils


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

    parser.add_argument('-F', '--format',
                        default='yaml', choices=['yaml', 'json'],
                        help='Workflow/Plugin file format')

    parser.add_argument('-s', '--strategy',
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

    strategy = getattr(launchers, args.strategy)

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

    if args.file:
        # Generate a new workflow
        wf = workflow.Workflow()
        for file_path in args.file:
            wf.add_component_from_file(file_path)

        # Load the project
        p = Project()

        # Run the workflow
        p.run(wf, strategy)

    elif args.plugin:
        # Generate a new workflow
        wf = workflow.Workflow()
        for plugin_id in args.plugin:
            wf.add_component(plugin_id)

        # Load the project
        p = Project()

        # Run the workflow
        p.run(wf, strategy)

    elif args.workflow:
        # Load the workflow (built beforehand and saved to a file)
        data = utils.load_struct(args.path, args.format)
        wf = workflow.from_node_link_data(data)

        # Load the project
        p = Project()

        # Run the workflow
        p.run(wf, strategy=strategy)

    else:
        raise ValueError(opts_err)
