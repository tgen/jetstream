"""Run a Jetstream workflow.

Template variable arguments should follow the syntax: ``--<key> <value>``.
The key must start with two hyphens and the value is the following argument. The
variable type can be explicitly set with the syntax ``--<type>:<key> <value>``.
Variables with no type declared will be loaded as strings.

If the variable type is "file" the value will be passed to
``jetstream.data_loaders``, handled to the extension. All other types will
evaluated by the appropriate type function.

"""
import sys
import logging
import argparse
import jetstream
from jetstream.cli import shared
from jetstream.backends import LocalBackend, SlurmBackend

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream pipelines',
        description=__doc__.replace('``', '"'),
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('templates', help='Template path', nargs='+')

    parser.add_argument('-t', '--search-path', action='append', default=None,
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

    parser.add_argument('--autosave', default=0, type=int,
                        help='Automatically save the workflow during run.')

    parser.add_argument('-w', '--workflow',
                        help='Run an existing workflow. Generally used for retry/resume '
                             'when a workflow run fails.')

    parser.add_argument('--method', choices=['retry', 'resume'], default='retry',
                        help='Method to use when restarting existing workflows')

    parser.add_argument('--retry', dest='method', action='store_const', const='retry',
                        help='Reset "failed" and "pending" tasks before starting')

    parser.add_argument('--resume',  dest='method', action='store_const', const='resume',
                        help='Reset "pending" tasks before starting')

    parser.add_argument('--max-forks', default=None, type=int,
                        help='Override the default fork limits of the task '
                             'backend.')

    return parser


def main(args=None):
    parser = arg_parser()
    args, remaining = parser.parse_known_args(args)
    log.debug(args)

    if args.workflow:
        workflow = jetstream.load_workflow(args.workflow)
        if args.method == 'retry':
            workflow.retry()
        else:
            workflow.resume()
    else:
        data = vars(shared.parse_kvargs(
            args=remaining,
            type_separator=args.kvarg_separator
        ))

        log.debug('Template render data: {}'.format(data))

        templates = jetstream.render_templates(
            *args.templates,
            data=data,
            search_path=args.search_path
        )

        if args.render_only:
            for t in templates:
                print(t)
            sys.exit(0)

        workflow = jetstream.build_workflow('\n'.join(templates))

    log.debug('Workflow data: {}'.format(workflow))

    if args.build_only:
        print(workflow.to_yaml())
        sys.exit(0)

    if args.backend == 'slurm':
        backend = SlurmBackend()
    else:
        backend = LocalBackend()

    runner = jetstream.AsyncRunner(
        backend=backend,
        max_forks=args.max_forks,
        autosave=args.autosave
    )

    rc = runner.start(workflow=workflow)
    runner.close()

    sys.exit(rc)


if __name__ == '__main__':
    main()
