"""Command line utility for managing plugins """
import sys
import argparse
import logging
from jetstream import plugins, utils

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.set_defaults(action=main)

    parser.add_argument('action',
                        choices=['ls', 'list', 'update', 'clone', 'rm',
                                 'remove', 'get'])

    parser.add_argument('-s', '--script-only',
                        action='store_true', default=False)

    parser.add_argument('-f', '--format',
                        default='yaml', help=argparse.SUPPRESS)

    parser.add_argument('plugin_id', nargs='?')

    return parser


def ls(plugin_id=None):
    if plugin_id is not None:
        revs = plugins.list_revisions(plugin_id)
        print(utils.yaml_dumps(revs))
    else:
        print(utils.yaml_dumps(plugins.ls()))


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    if args.action in ('update',):
        plugins.update()

    elif args.action in ('clone',):
        plugins.clone(args.plugin_id)

    elif args.action in ('rm', 'remove'):
        plugins.remove(args.plugin_id)

    elif args.action in ('get',):
        if args.plugin_id:
            plugin = plugins.get_plugin(
                args.plugin_id,
                script_only=args.script_only
            )
            if args.script_only:
                print(plugin)
            else:
                text = utils.yaml_dumps(plugin)
                print(text)

        else:
            print('No plugin id given, available plugins: ')
            ls(args.plugin_id)
            sys.exit(1)
    else:
        ls(args.plugin_id)
