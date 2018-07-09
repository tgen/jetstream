"""Directly render templates from file/string.

This command inherits most of its behavior from jetstream_pipelines, but
loads templates directly from a path or string instead of searching for them
in a template directory.

Note: This command is primarily intended for debugging template issues. To
render, build, and run a workflow in one step use "jetstream_pipelines".
"""
import logging
import jetstream
from jetstream.cli.subcommands.pipelines import arg_parser as base_parser
from jetstream.cli import kvargs

log = logging.getLogger(__name__)


def arg_parser():
    parser = base_parser()
    parser.epilog = __doc__

    parser.add_argument('-s', '--string', action='store_true', default=False,
                        help='Render a templete from a string instead of '
                             'reading a file.')

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    try:
        project = jetstream.Project()
    except Exception as e:
        log.exception(e)
        project = dict()

    if args.string:
        template_text = args.template
    else:
        with open(args.template, 'r') as fp:
            template_text = fp.read()

        log.critical("Loaded template: {}".format(args.template))

    template = jetstream.env.from_string(template_text)

    additional_data = kvargs.parse(args.kvargs)

    print(template.render(project=project, **vars(additional_data)))
