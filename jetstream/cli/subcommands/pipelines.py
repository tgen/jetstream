"""Run a Jetstream pipeline

All arguments following the name of the pipeline are considered template
variable data and will be ignored by the argument parser initially. Any command
options listed below *must* be given *before* the template name.
However, logging arguments (see ``jetstream -h``) are special and can be
included anywhere.

*Template variable data:*

Template variable data is usually saved to files in ``<project>/config``, but
command arguments can also be used to pass variable data to templates.
Variables are specified as argument pairs: ``--<key> <value>``. The key must
start with two hyphens and the value is the next argument.

Variables can also be typed with the syntax ``--<type>:<key> <value>``.
Some popular types are "file", "json", and "yaml". Files will be handled with
``jetstream.data_loaders`` according to their extension. All others will
evaluated by the appropriate loader function. Variables with no type declared
will be loaded as strings.

"""
import os
import sys
import logging
import argparse
import jetstream
from jetstream.cli import kvargs

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream pipelines',
        description=__doc__.replace('``', '"'),
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('template', help='Template name')

    # parser.add_argument('kvargs', nargs=argparse.REMAINDER,
    #                     help='Arguments following the workflow name are parsed '
    #                          'as arbitrary "--key value" pairs (see help)')

    parser.add_argument('--kvarg-separator', default=':',
                        help='Specify an alternate separator for kvargs')

    parser.add_argument('-r', '--render-only', action='store_true')

    parser.add_argument('--backend', default='LocalBackend',
                        help=argparse.SUPPRESS)

    return parser


def build_workflow(run, tasks):
    workflow = jetstream.workflows.build_workflow(tasks)
    workflow.save(os.path.join(run.path, 'workflow'))
    return workflow


def render_tasks(run, template, project, additional_data):
    tasks = template.render(project=project, **vars(additional_data))
    with open(os.path.join(run.path, 'tasks'), 'w') as fp:
        fp.write(tasks)
    return tasks


def get_template(run, template_name):
    template = jetstream.env.get_template_with_source(template_name)
    with open(os.path.join(run.path, 'template'), 'w') as fp:
        jetstream.utils.yaml.dump(template.source, stream=fp)
    return template


def main(args=None):
    parser = arg_parser()
    args, unknown = parser.parse_known_args(args)
    log.debug(args)

    kvargs_data = kvargs.parse(
        args=unknown,
        type_separator=args.kvarg_separator
    )

    p = jetstream.Project()

    if args.render_only:
        text = p.render(
            template=args.template,
            additional_data=vars(kvargs_data)
        )

        print(text)

    else:
        backend = getattr(jetstream, args.backend)

        rc = p.run(
            template=args.template,
            additional_data=vars(kvargs_data),
            backend=backend
        )

        sys.exit(rc)


if __name__ == '__main__':
    main()
