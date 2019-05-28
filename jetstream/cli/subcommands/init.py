"""Create or reinitialize a project

This command is used to create a new Jetstream project directory. If no
path is given, the current directory will be initialized. If config data
options are given (-c/--config/--config-file), they will be added to the
project config file.
"""
import os
import logging
import jetstream
from jetstream.cli import add_config_options_to_parser

log = logging.getLogger(__name__)

def arg_parser(parser):
    parser.add_argument(
        'path',
        nargs='?',
        default=os.getcwd(),
        help='Path to a initialize a project'
    )

    parser.add_argument(
        '-f', '--force',
        action='store_true',
        help='Force overwrite of project.yaml'
    )

    parser.add_argument(
        '--project-id',
        help='Force a project ID instead of using letting it be generated '
             'automatically'
    )

    add_config_options_to_parser(parser)
    return parser


def main(args):
    log.debug(f'{__name__} {args}')
    p = jetstream.Project(args.path)

    if p.exists():
        if args.force:
            log.info('Project already exists - force reinitializing')
            p.init(id=args.project_id, config=args.config)
        else:
            log.info('Project already exists - updating config')
            p.save_config(args.get_config)
    else:
        log.info('Initializing project')
        p.init(id=args.project_id, config=args.config)

    log.info(p)
    return p
