"""Command line utility for launching a single plugin component """
import argparse
import json
import logging
from jetstream import launchers, plugins, workflow, Project
from sys import exit

log = logging.getLogger(__name__)

def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('path',
                        help='Path to a plugin file, workflow, or plugin ID '
                             'in the format: "plugin/path<:revision>')

    parser.add_argument('-f', '--file',
                        action='store_true', default=False,
                        help='Launch plugin from file instead of repo')

    parser.add_argument('-w', '--workflow',
                        action='store_true', default=False,
                        help='Launch a workflow')

    parser.add_argument('-s', '--strategy',
                        help=argparse.SUPPRESS,
                        default='default',
                        choices=['dry', 'default'])

    parser.add_argument('--no-records',
                        action='store_true', default=False)

    return parser


def launch_single_plugin(plugin, strategy, no_records=False):
    #TODO HOW TO LAUNCH A PLUGIN FROM A FILE AND ADD TO WORKFLOW?
    if no_records:
        result = strategy(plugin)
        print(json.dumps(result, indent=4))


def main(args):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    strategy = getattr(launchers, args.strategy)

    if args.file and args.workflow:
        parser.print_help()
        print("'-f/--file' and '-w/--workflow' are mutually exclusive")
        exit(1)

    if args.file:
        plugin = plugins.load(args.path)
    elif args.workflow:
        # Load the workflow (built beforehand and saved to a file)
        with open(args.path, 'r') as fp:
            data = json.load(fp)

        wf = workflow.from_json(data)

        # Start the runner (there may be more than one of these if performance
        # becomes a concern)
        p = Project()
        p.run(wf, strategy=strategy)
    else:
        plugin = plugins.get_plugin(args.path)


