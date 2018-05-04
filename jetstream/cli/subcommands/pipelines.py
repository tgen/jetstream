"""Run a Jetstream pipeline

All arguments following the name of the pipeline are considered template data
values and will be ignored by the argument parser initially. Any options listed
below must be given BEFORE the workflow name (first positional argument) The
key must start with two hyphens and the value is the next argument. Values can
be JSON strings.

    jetstream pipelines <options> <template> [ --<key> <value> ]

"""
import sys
import logging
import argparse
import jetstream

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream pipelines',
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

    parser.add_argument('--kvarg-separator', default=':', help=argparse.SUPPRESS)

    return parser


# Loader functions for typed arbitrary arguments
argtype_fns = {
    "default": str,
    "file": jetstream.project.load_data_file,
    "json": jetstream.utils.json_loads,
    "yaml": jetstream.utils.yaml_loads,
}


def reparse_aribitrary(args, type_separator=':'):
    """Reparses sequence of arbitrary arguments "--<type>:<key> <value>"

    This works by first building an argument parser specifically for the
    arguments present in the list. First we look for any items that
    start with '--', then adding an argument to the parser with the given
    key and type (type is optional, str is the default).

    After building the parser, the args are parsed and namespace is returned
    as a dictionary. """
    parser = argparse.ArgumentParser(add_help=False)

    for arg in args:
        if arg.startswith('--'):

            if type_separator in arg:
                argtype, _, key = arg.lstrip('-').partition(type_separator)
            else:
                argtype = 'default'
                key = arg.lstrip('-')

            log.debug('Adding parser key: {} type:{}'.format(key, argtype))
            fn = argtype_fns[argtype]
            parser.add_argument(arg, type=fn, dest=key)

    namespace = parser.parse_args(args)
    return vars(namespace)


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.critical(args)

    env = jetstream.template_env(
        strict=args.strict,
        include_project_templates=args.project_templates,
    )

    project = jetstream.Project()
    t = env.get_template(args.template)
    log.critical("Loaded template: {}".format(t.filename))

    # Any arguments following the workflow template name are parsed with
    # json_allowed and available as variables when rendering the template.
    # "project" is reserved as the namespace for the project data.
    kwargs = reparse_aribitrary(
        args=args.kvargs,
        type_separator=args.kvarg_separator
    )
    log.debug('Project config data:\n{}'.format(project.config))
    log.debug('Other template data:\n{}'.format(kwargs))
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
