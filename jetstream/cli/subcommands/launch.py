"""Command line utility for launching a single plugin component """
import argparse
import json
import logging
from jetstream import launch, plugins

log = logging.getLogger(__name__)

def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('plugin_id',
                        help='Plugin identifier in the format: "plugin/path<:revision>')

    parser.add_argument('--file',
                        action='store_true', default=False,
                        help='Load plugin from file instead of repo')

    parser.add_argument('-s', '--strategy',
                        default='default',
                        choices=['dry', 'default'])
    return parser


def main(args):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    launcher = getattr(launch, args.strategy)

    if args.file:
        plugin = plugins.load(args.plugin_id)
    else:
        plugin = plugins.get_plugin(args.plugin_id)

    result = launcher(plugin)
    print(json.dumps(result, indent=4))
