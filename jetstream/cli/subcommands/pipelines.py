"""Run a Jetstream pipeline

Template variable data is usually saved as files in ``<project>/config``, but 
command arguments can also be used to pass variable data to templates. All 
arguments remaining after parsing the command line (arguments that are not 
listed in the help section) will be treated as template variable data (kvargs):

Template variable arguments should follow the syntax: ``--<key> <value>``. 
The key must start with two hyphens and the value is the following argument. The 
variable type can be explicitly set with the syntax ``--<type>:<key> <value>``.
Variables with no type declared will be loaded as strings.

If the variable type is "file" the value will be ``jetstream.data_loaders``, which
handles files according to their extension. All other types will evaluated by the 
appropriate type function. 

"""
import os
import sys
import logging
import argparse
import subprocess
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

    parser.add_argument('--kvarg-separator', default=':',
                        help='Specify an alternate separator for kvargs')

    parser.add_argument('-r', '--render-only', action='store_true',
                        help='Just render the template and print to stdout')

    parser.add_argument('--backend', choices=['local', 'slurm'], default='local',
                        help='Specify the runner backend (default: local)')

    parser.add_argument('--logging-interval', default=300, type=int,
                        help='Time between workflow status updates')

    parser.add_argument('--max-concurrency', default=None, type=int,
                        help='Override the concurrency limits of the task backend.')
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
        type_separator=args.kvarg_separator)

    p = jetstream.Project()

    if args.render_only:
        text = p.render(
            template=args.template,
            additional_data=vars(kvargs_data))

        print(text)

    else:
        if args.backend == 'local':
            backend = jetstream.LocalBackend(max_subprocess=args.max_concurrency)  
        elif args.backend == 'slurm':
            backend = jetstream.SlurmBackend(max_jobs=args.max_concurrency)

        backend_class = getattr(jetstream, args.backend)
        backend = backend_class()

        rc = p.run(
            template=args.template,
            additional_data=vars(kvargs_data),
            backend=backend,
            logging_interval=args.logging_interval)

        sys.exit(rc)


if __name__ == '__main__':
    main()
