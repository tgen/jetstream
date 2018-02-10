"""Command line utility for managing plugins """
import json
from jetstream import plugins

def arg_parser(subparser):
    parser = subparser.add_parser('plugins', description=__doc__)
    parser.set_defaults(action=main)

    parser.add_argument('subaction',
                        choices=['ls', 'list', 'update', 'add', 'rm', 'get',
                                 'get-script'])

    parser.add_argument('plugin_id', nargs='?')


def main(args):
    if args.subaction in ('update',):
        plugins.sync()
    elif args.subaction in ('add',):
        plugins.clone(args.plugin_id)
    elif args.subaction in ('rm',):
        plugins.remove(args.plugin_id)
    elif args.subaction in ('get',):
        plugin_data = plugins.get_plugin(args.plugin_id)
        print(json.dumps(plugin_data, indent=4))
    elif args.subaction in ('get-script',):
        plugin_data = plugins.get_plugin(args.plugin_id)
        print(plugin_data['script'])
    else:
        # ls, list are default action, fallback to this
        if args.plugin_id:
            print(plugins.list_revisions(args.plugin_id))
        else:
            print(plugins.ls())
