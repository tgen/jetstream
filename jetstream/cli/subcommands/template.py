"""Build a Jetstream workflow from template files.

Template variable arguments should follow the syntax: ``--<key> <value>``.
The key must start with two hyphens and the value is the following argument. The
variable type can be explicitly set with the syntax ``--<type>:<key> <value>``.
Variables with no type declared will be loaded as strings.

If the variable type is "file" the value will be passed to
``jetstream.data_loaders``, All other types will evaluated by their type
function.

"""
import sys
import logging
import argparse
import jetstream
from jetstream import utils
from jetstream.cli import shared

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream build',
        description=__doc__.replace('``', '"'),
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('template', help='Template path')

    parser.add_argument('out', help='Where to save the workflow')

    parser.add_argument('-c', '--config',
                        help='Load a single config file instead of loading '
                             'project files')

    parser.add_argument('-r', '--render-only', action='store_true',
                        help='Just render the template and print to stdout')

    parser.add_argument('-t', '--search-path', action='append', default=None,
                        help='Manually configure the template search. This '
                             'argument can be used multiple times.')

    parser.add_argument('--format', choices=[None, 'yaml', 'json', 'pickle'],
                        default=None)

    parser.add_argument('--kvarg-separator',
                        help='Specify an alternate separator for kvargs')

    return parser


def build_context(config, kvargs, kvarg_separator=':'):
    try:
        project = jetstream.Project()
    except jetstream.projects.NotAProject:
        project = None

    if project and config:
        context = utils.yaml_load(config)
        context['project'] = project
    elif project:
        context = project.config
        context['project'] = project
    elif config:
        context = utils.yaml_load(config)
        context['project'] = None
    else:
        context = {}

    kvarg_data = vars(shared.parse_kvargs(
        args=kvargs,
        separator=kvarg_separator
    ))

    context.update(kvarg_data)
    return context


def main(args=None):
    parser = arg_parser()
    args, remaining = parser.parse_known_args(args)
    log.debug(args)

    context = build_context(
        config=args.config,
        kvargs=remaining,
        kvarg_separator=args.kvarg_separator
    )

    log.debug('Template render context: {}'.format(context))

    render = jetstream.render_template(
        args.template,
        data=context,
        search_path=args.search_path
    )

    if args.render_only:
        print(render)
    else:
        workflow = jetstream.build_workflow(render)
        log.debug('Workflow data: {}'.format(workflow))

        workflow.save(args.out, format=args.format)


if __name__ == '__main__':
    main()
