"""Command line utility for launching a single plugin component """
import json
from jetstream import launch, plugins


def arg_parser(subparser):
    parser = subparser.add_parser('launch', description=__doc__)
    parser.set_defaults(action=main)

    parser.add_argument('plugin_id',
                        help='Plugin identifier in the format: "plugin/path<:revision>')

    parser.add_argument('--file',
                        action='store_true', default=False,
                        help='Load plugin from file instead of repo')

    parser.add_argument('-s', '--strategy',
                        default='default',
                        choices=['dry', 'default'])


def main(args):
    launcher = getattr(launch, args.strategy)

    if args.file:
        plugin = plugins.load(args.file)
    else:
        plugin = plugins.get_plugin(args.plugin_id)

    result = launcher(plugin)
    print(json.dumps(result, indent=4))
