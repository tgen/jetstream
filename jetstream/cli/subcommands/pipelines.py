"""Run a Jetstream workflow inside a project


Template variable data is usually saved as files in ``<project>/config``, but 
command arguments can also be used to pass variable data to templates. All 
arguments remaining after parsing the command line (arguments that are not 
listed in the help section) will be treated as template variable data (kvargs):

Template variable arguments should follow the syntax: ``--<key> <value>``. 
The key must start with two hyphens and the value is the following argument. The 
variable type can be explicitly set with the syntax ``--<type>:<key> <value>``.
Variables with no type declared will be loaded as strings.

If the variable type is "file" the value will be ``jetstream.data_loaders``,
which handles files according to their extension. All other types will
evaluated by the appropriate type function.

"""
import sys
import logging
import argparse
import jetstream
from jetstream.cli import shared

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream pipelines',
        description=__doc__.replace('``', '"'),
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('template', help='Template path', nargs='+')

    parser.add_argument('-t', '--template-search-path', action='append',
                        help='Manually configure the template search. This '
                             'argument can be used multiple times.')

    parser.add_argument('--kvarg-separator', default=':',
                        help='Specify an alternate separator for kvargs')

    parser.add_argument('-r', '--render-only', action='store_true',
                        help='Just render the template and print to stdout')

    parser.add_argument('-b', '--build-only', action='store_true',
                        help='Just build the workflow and print to stdout')

    parser.add_argument('--backend', choices=['local', 'slurm'], default='local',
                        help='Specify the runner backend (default: local)')

    parser.add_argument('--logging-interval', default=60, type=int,
                        help='Time between workflow status updates')

    parser.add_argument('--max-forks', default=None, type=int,
                        help='Override the fork limits of the task backend.')

    return parser


def main(args=None):
    parser = arg_parser()
    args, unknown = parser.parse_known_args(args)
    log.debug(args)

    kvargs_data = shared.parse_kvargs(
        args=unknown,
        type_separator=args.kvarg_separator)

    p = jetstream.Project()

    tasks = list()

    for path in args.template:
        template = shared.load_template(path, args.template_search_path)
        rendered = template.render(project=p, **vars(kvargs_data))

        if args.render_only:
            print(rendered)
            continue

        loaded = jetstream.utils.yaml_loads(rendered)
        tasks.extend(loaded)

    if args.render_only:
        pass
    elif args.build_only:
        wf = jetstream.workflows.build_workflow(tasks)
        print(wf.pretty())
    else:
        wf = jetstream.workflows.build_workflow(tasks)

        if args.backend == 'slurm':
            backend = jetstream.SlurmBackend(max_forks=args.max_forks)
        else:
            backend = jetstream.LocalBackend(max_forks=args.max_forks)

        rc = p.run(wf, backend=backend, logging_interval=args.logging_interval)
        sys.exit(rc)


if __name__ == '__main__':
    main()
