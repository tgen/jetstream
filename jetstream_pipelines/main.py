"""Run a Jetstream pipeline

All arguments following the name of the pipeline are considered template data
values and will be ignored by the argument parser initially. Any options listed
below must be given BEFORE the workflow name (first positional argument) The
key must start with two hyphens and the value is the next argument. Values can
be JSON strings.

    jetstream_pipelines <options> <template> [ --<key> <value> ]

"""
import sys
import json
import logging
import argparse
import jetstream
import jetstream_pipelines

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('template', help='Template name')

    parser.add_argument('-b', '--build-only', action='store_true',
                        help='Just build the workflow and print to stdout.')

    parser.add_argument('-r', '--render-only', action='store_true',
                        help='Just render the template and print to stdout.')

    parser.add_argument('--ignore-undefined', dest='strict',
                        action='store_false',
                        help='Suppress errors normally raised when a workflow '
                             'variable is undefined')

    parser.add_argument('--no-project-templates', dest='project_templates',
                        action='store_false',
                        help='Ignore templates in the current project')

    parser.add_argument('kvargs', nargs=argparse.REMAINDER,
                        help='Arguments following the workflow name are parsed '
                             'as arbitrary "--key value" pairs (see help)')

    return parser



def reparse_aribitrary(remainder):
    parser = argparse.ArgumentParser(add_help=False)
    for arg in remainder:
        # TODO if the value starts with a dash, what do?
        if arg.startswith(("-", "--")):
            parser.add_argument(arg, type=json_allowed)
    return vars(parser.parse_args(remainder))


def json_allowed(value):
    """Given any value this will attempt to parse as json and return the
    resulting object. If parsing fails, it will fall back to a string."""
    try:
        return json.loads(value)
    except json.decoder.JSONDecodeError:
        log.debug('JSON parsing failed, treating as string: {}'.format(value))
        return str(value)


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.critical(args)

    env = jetstream_pipelines.env(
        strict=args.strict,
        include_project_templates=args.project_templates,
    )

    project = jetstream.Project()
    t = env.get_template(args.template)
    log.critical("Loaded template: {}".format(t.filename))

    # Any arguments following the workflow template name are parsed with
    # json_allowed and available as variables when rendering the template.
    # "project" is reserved as the namespace for the project data.
    kwargs = reparse_aribitrary(args.kvargs)
    rendered_template = t.render(project=project, **kwargs)

    if args.render_only:
        print(rendered_template)
        sys.exit(0)

    # Node properties are described in jetstream.workflows.spec
    nodes = jetstream.utils.yaml_loads(rendered_template)
    wf = jetstream.workflows.build_workflow(nodes)

    if args.build_only:
        print(wf)
        sys.exit(0)

    jetstream.workflows.run_workflow(wf)

    failures = [s for s in wf.status() if s == 'failed']
    if failures:
        log.critical('Error: Some tasks failed! {}'.format(failures))
        sys.exit(1)

if __name__ == '__main__':
    main()
