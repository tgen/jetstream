"""
"""
import os
import logging
import argparse
import jetstream
from jetstream.cli import shared

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream pipelines',
        description=__doc__.replace('``', '"'),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        'name',
        help='Pipeline name'
    )

    parser.add_argument(
        '--project',
        help='Path to a project'
    )

    parser.add_argument(
        '--separator',
        default=':',
        help='Set an alternate separator for kvargs variables'
    )

    parser.add_argument(
        '--pipelines-dir',
        default=os.environ.get('JETSTREAM_PIPELINES'),
        help='Override path to the pipelines dir',
    )

    runner = parser.add_argument_group(
        title='Runner Configuration'
    )

    runner.add_argument(
        '--save-interval',
        type=int,
        default=30,
        help='Frequency, in seconds, that the workflow file will be saved to '
             'disk. (Default: 3600)')

    runner.add_argument(
        '--backend',
        choices=['local', 'slurm', None],
        default=None,
        help='Specify the runner backend (Default: $JETSTREAM_BACKEND or local)'
    )

    runner.add_argument(
        '--local',
        dest='backend',
        action='store_const',
        const='local'
    )

    runner.add_argument(
        '--slurm',
        dest='backend',
        action='store_const',
        const='slurm'
    )

    runner.add_argument(
        '--run-id',
        help='Give this run a specific ID instead of randomly generating one.'
    )

    runner.add_argument(
        '--max-forks',
        default=None,
        type=int,
        help='Override the default fork limits of the runner. By default, '
             'this will be set to a conservative fraction of the total thread '
             'limit on runner host.')

    runner.add_argument(
        '--method',
        choices=['retry', 'resume', 'reset'],
        default='retry',
        help='Method to use when running existing workflows. This parameter '
             'will determine which tasks are reset prior to starting the run.'
    )

    return parser


def main(args=None):
    parser = arg_parser()
    args, remaining = parser.parse_known_args(args)
    log.debug(args)

    pipeline, constants = jetstream.pipelines.get(args.name)

    if constants is not None:
        remaining.extend(['--file:constants', str(constants.resolve())])

    # If the pipeline has a bin directory, prepend the env PATH
    if pipeline.bin:
        bin_path = pipeline.path.joinpath(pipeline.manifest['bin'])
        os.environ['PIPELINE_BIN'] = str(bin_path)
        os.environ['PATH'] = f'{bin_path}:{os.environ["PATH"]}'

    context = shared.load_context(
        project=args.project,
        kvargs=remaining,
        kvarg_separator=args.separator
    )

    workflow = jetstream.render_template(
        path=pipeline.main,
        context=context
    )

    if args.project:
        project = jetstream.Project(path=args.project)
        existing_wf = project.workflow()

        if existing_wf is not None:
            workflow = jetstream.workflows.mash(existing_wf, workflow)

    if args.method == 'retry':
        workflow.retry()
    elif args.method == 'resume':
        workflow.resume()
    elif args.method == 'reset':
        workflow.reset()

    runner = jetstream.Runner(
        backend=args.backend,
        max_forks=args.max_forks,
        autosave=args.save_interval,
    )

    runner.start(
        workflow=workflow,
        project=args.project,
        run_id=args.run_id,
    )

if __name__ == '__main__':
    main()
