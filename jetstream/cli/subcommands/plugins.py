"""Command line utility for managing plugins """
import sys
import argparse
import logging
from jetstream import plugins

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.set_defaults(action=main)

    parser.add_argument('action',
                        choices=['ls', 'list', 'update', 'clone', 'rm',
                                 'remove', 'get', 'get-script'])

    parser.add_argument('plugin_id', nargs='?')

    return parser


def ls(plugin_id=None):
    if plugin_id is not None:
        print(plugins.list_revisions(plugin_id))
    else:
        [print(p) for p in plugins.ls()]


def main(args):
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
            print(plugins.get_plugin(args.plugin_id, raw=True))
        else:
            print('No plugin id given, available plugins: ')
            ls(args.plugin_id)
            sys.exit(1)

    elif args.action in ('get-script',):
        if args.plugin_id:
            plugin_data = plugins.get_plugin(args.plugin_id)
            print(plugin_data['script'])
        else:
            print('No plugin id given, available plugins: ')
            ls(args.plugin_id)
            sys.exit(1)

    else:
        ls(args.plugin_id)
