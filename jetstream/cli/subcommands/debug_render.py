"""Directly render templates from file/string.

This command inherits most of its behavior from jetstream_pipelines, but
loads templates directly from a path or string instead of searching for them
in a template directory.

Note: This command is primarily intended for debugging template issues. To
render, build, and run a workflow in one step use "jetstream_pipelines".
"""
import logging
import jetstream
import jetstream_pipelines
from jetstream_pipelines.main import arg_parser, reparse_aribitrary

log = logging.getLogger(__name__)


def main(args=None):
    parser = arg_parser()
    parser.epilog = __doc__
    parser.add_argument('-s', '--string', action='store_true', default=False,
                        help='Render a templete from a string instead of '
                             'reading a file.')

    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    project = jetstream.Project()

    # Jinja manages templates, we just need to add any additional directories
    # to the search path, and decide if we want strict rendering.
    env = jetstream_pipelines.env(strict=args.strict)

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
