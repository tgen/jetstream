"""Directly render templates from file/string.

This command inherits most of its behavior from jetstream_pipelines, but
loads templates directly from a path or string instead of searching for them
in a template directory.

Note: This command is primarily intended for debugging template issues. To
render, build, and run a workflow in one step use "jetstream_pipelines".
"""
import logging

import jetstream
from jetstream.cli.subcommands.pipelines import arg_parser, reparse_aribitrary

log = logging.getLogger(__name__)


def main(args=None):
    parser = arg_parser()
    parser.epilog = __doc__
    parser.add_argument('-s', '--string', action='store_true', default=False,
                        help='Render a templete from a string instead of '
                             'reading a file.')

    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    try:
        project = jetstream.Project()
    except Exception as e:
        log.exception(e)
        project = dict()

    env = jetstream.template_env(
        strict=args.strict,
        include_project_templates=args.project_templates
    )

    if args.string:
        template_text = args.template
    else:
        with open(args.template, 'r') as fp:
            template_text = fp.read()
        log.critical("Loaded template: {}".format(args.template))

    template = env.from_string(template_text)

    # Any arguments following the workflow template name are parsed with
    # json_allowed and available as variables when rendering the template.
    # "project" is reserved as the namespace for the project data.
    kwargs = reparse_aribitrary(args.kvargs)
    rendered_template = template.render(project=project, **kwargs)

    print(rendered_template)
