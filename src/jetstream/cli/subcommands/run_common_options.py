import jetstream
from jetstream.cli import add_config_options_to_parser


def add_arguments(parser):
    group = parser.add_argument_group(
        'execution options',
        description='These options are available for any of the commands that build or '
                    'execute workflows.'
    )

    group.add_argument(
        '-o', '--out',
        help='path to save the render, build, or running workflow progress '
             '[automatically set if working with a project]'
    )

    group.add_argument(
        '-m', '--mash-only',
        action='store_true',
        help='just render the template, build the workflow, mash with an existing workflow, and stop'
    )

    group.add_argument(
        '-b', '--build-only',
        action='store_true',
        help='just render the template, build the workflow, and stop'
    )

    group.add_argument(
        '-r', '--render-only',
        action='store_true',
        help='just render the template and stop'
    )

    group.add_argument(
        '--backend',
        choices=jetstream.settings['backends'].flatten(),
        default=jetstream.settings['backend'].get(str),
        help='runner backend name used for executing tasks [%(default)s]'
    )

    group.add_argument(
        '--reset-method',
        choices=['retry', 'resume', 'reset'],
        default='retry',
        help='controls which tasks are reset prior to starting the run - '
             '"retry": pending and failed, "resume": pending, or "reset": '
             'all [%(default)s]'
    )

    group.add_argument(
        '-w', '--existing-workflow',
        help='path to an existing workflow file that will be merged into run '
             '[automatically set if working with a project]'
    )

    group.add_argument(
        '-s', '--search-path',
        action='append',
        default=[],
        help='directory to add to search path for loading templates and pipelines,'
             'this can be used multiple times'
    )

    add_config_options_to_parser(parser)
