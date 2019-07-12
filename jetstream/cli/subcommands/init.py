"""Create or update a project.

This command is used to create a new Jetstream project directory. If no
path is given, the current directory will be initialized. If config variables
arguments are given (-c/--config -C/--config-file), they will be added to the
project manifest file.
"""
import os
import logging
import jetstream
from jetstream.cli import add_config_options_to_parser

log = logging.getLogger('jetstream.cli')

def arg_parser(parser):
    parser.add_argument(
        'path',
        nargs='?',
        default=os.getcwd(),
        help='directory where project files should be created/updated'
    )

    parser.add_argument(
        '-f', '--force',
        action='store_true',
        help='force overwrite of project.yaml'
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

    paths = jetstream.projects.ProjectPaths(args.path)
    if paths.exists():
        if args.force:
            log.info('Project already exists - force reinitializing')
            p = jetstream.projects.init(
                args.path,
                id=args.project_id,
                config=args.config
            )
        else:
            log.info('Project already exists - updating config')
            p = jetstream.Project(args.path)
            p.update_index(args.config)
    else:
        log.info('Initializing project')
        p = jetstream.projects.init(
            args.path,
            id=args.project_id,
            config=args.config
        )

    log.info(p)
