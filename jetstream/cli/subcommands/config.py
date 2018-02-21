""" This module contains the cli interface code for the config file utility.
"""


def arg_parser(subparser):
    parser = subparser.add_parser('config', description=__doc__)
    parser.set_defaults(action=main)

    parser.add_argument('subaction',
                        choices=['legacy-convert', 'validate'])

    parser.add_argument('--format', default='json',
                        choices=('json', 'yaml'))
    parser.add_argument('path')


def main(args):
    if args.subaction in ('legacy-convert',):
        from jetstream import config
        c = config.legacy.load(args.path)
        print(config.serialize(c, format=args.format))
    elif args.subaction in ('validate',):
        from jetstream import config
        c = config.load(args.path, format=args.format)
        config.validate(c)
