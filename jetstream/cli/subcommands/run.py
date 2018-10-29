"""Run a Jetstream workflow.

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
from jetstream.cli import shared
from jetstream.backends import LocalBackend, SlurmBackend

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream pipelines',
        description=__doc__.replace('``', '"'),
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('workflow', help='Workflow path')

    parser.add_argument('--autosave', type=int,
                        default=jetstream.settings['autosave'],
                        help='Automatically save the workflow during run.')

    parser.add_argument('--backend', choices=['local', 'slurm'],
                        default=jetstream.settings['backend'],
                        help='Specify the runner backend (default: local)')

    parser.add_argument('--format', choices=[None, 'yaml', 'json', 'pickle'],
                        default=None)

    parser.add_argument('--method', choices=['retry', 'resume', 'reset'],
                        default='retry',
                        help='Method to use when restarting workflows. This '
                             'will determine which tasks are reset prior to '
                             'starting.')

    parser.add_argument('--resume', dest='method', action='store_const',
                        const='resume',
                        help='Reset "pending" tasks before starting')

    parser.add_argument('--retry', dest='method', action='store_const',
                        const='retry',
                        help='Reset "failed" and "pending" tasks before starting')

    parser.add_argument('--run-id',
                        help='Give this run a specific ID instead of randomly '
                             'generating one.')

    parser.add_argument('--max-forks', default=None, type=int,
                        help='Override the default fork limits of the task '
                             'backend.')

    parser.add_argument('--check', action='store_true', dest='check_workflow',
                        default=True,
                        help='Check workflow for failed tasks after run.')

    parser.add_argument('--no-check', action='store_false', dest='check_workflow',
                        help='Ignore failed tasks after run. Runtime errors '
                             '(problems launching tasks, or runner errors) will '
                             'still cause a non-zero exit status.')

    return parser


def main(args=None):
    parser = arg_parser()
    args, remaining = parser.parse_known_args(args)
    log.debug(args)

    workflow = jetstream.load_workflow(args.workflow, args.format)

    if args.method == 'retry':
        workflow.retry()
    elif args.method == 'resume':
        workflow.resume()
    elif args.method == 'reset':
        workflow.reset()

    if args.backend == 'slurm':
        backend = SlurmBackend(max_concurrency=9002)
    else:
        backend = LocalBackend()

    runner = jetstream.Runner(
        backend=backend,
        max_concurrency=args.max_forks,
        autosave=args.autosave
    )

    runner.start(workflow=workflow, run_id=args.run_id)

    if args.check_workflow:
        shared.check_workflow_status(workflow)


if __name__ == '__main__':
    main()
