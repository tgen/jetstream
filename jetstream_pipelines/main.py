"""Run a Jetstream pipeline

All arguments following the name of the pipeline are considered template
data values and will be ignored by the argument parser initially. Any
options listed below must be given BEFORE the workflow name (first
positional argument)

Template Data "--key value" Arguments:

All arguments FOLLOWING the workflow name are parsed as key-value pairs,
and then passed to the template for rendering. The key must start with
two hyphens and the value is the next argument: "--key value". Values
can be JSON strings.

For example, if a workflow requires a variable "name":

    jetstream_pipelines bender --name bender

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

    parser.add_argument('name', help='Template name')

    parser.add_argument('-t', '--template-dir', action='append', default=[])

    parser.add_argument('--build-only', action='store_true',
                        help='Just build the workflow and print to stdout.')

    parser.add_argument('--render-only', action='store_true',
                        help='Just render the template and print to stdout.')

    parser.add_argument('--ignore-undefined', dest='strict',
                        action='store_false', default=True,
                        help='Suppress errors normally raised when a workflow '
                             'variable is undefined')

    parser.add_argument('--no-site-templates', dest='package_templates',
                        action='store_false', default=True,
                        help='Ignore templates installed in package data')

    parser.add_argument('--no-project-templates', dest='project_templates',
                        action='store_false', default=True,
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
    log.debug(args)

    project = jetstream.Project()
    template_name = args.name
    strict = args.strict
    additional_template_dirs = args.template_dir
    package_templates = args.package_templates
    project_templates = args.project_templates
    key_value_args = args.kvargs

    # TODO Template search path should be set as an environement variable
    # in order to propagate the settings to recursive workflows

    # Jinja manages templates, we just need to add any additional directories
    # to the search path, and decide if we want strict rendering.
    env = jetstream_pipelines.env(
        *additional_template_dirs,
        strict=strict,
        include_project_templates=project_templates,
        include_package_templates=package_templates
    )

    t = env.get_template(template_name)
    log.critical("Loaded template: {}".format(t.filename))

    # Any arguments following the workflow template name are parsed with
    # json_allowed and available as variables when rendering the template.
    # "project" is reserved as the namespace for the project data.
    kwargs = reparse_aribitrary(key_value_args)
    rendered_template = t.render(project=project, **kwargs)

    if args.render_only:
        print(rendered_template)
        sys.exit(0)

    # Node properties that are allowed in a workflow template are
    # described in jetstream.workflows.spec
    nodes = jetstream.utils.yaml_loads(rendered_template)
    wf = jetstream.workflows.build_workflow(nodes)

    if args.build_only:
        print(wf)
        sys.exit(0)

    jetstream.workflows.run_workflow(wf)


if __name__ == '__main__':
    main()
