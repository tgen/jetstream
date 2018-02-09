"""Command line utility for managing plugins """
from jetstream import reports


def arg_parser(subparser):
    parser = subparser.add_parser('report', description=__doc__)
    parser.set_defaults(action=main)

    parser.add_argument('project', nargs='*')


def main(args):
    print(args)