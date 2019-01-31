"""Create or reinitialize a project

This command can be used to create a new Jetstrema project with the
recommended project folders. If kvargs are given, they will be added to the
project config file."""
import os
import jetstream
from jetstream import log


def arg_parser(parser):
    parser.add_argument(
        'path',
        nargs='?',
        default=os.getcwd(),
        help='Path to a Jetstream project'
    )

    return parser


def main(args):
    log.debug(f'{__name__} {args}')
    p = jetstream.new_project(args.path, config=args.kvargs)
    # TODO append vs overwrite config and project file for existing projects?
    log.info(f'Initialized {p}')
