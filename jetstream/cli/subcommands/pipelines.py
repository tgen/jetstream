"""

"""
import os
import logging
import argparse
import jetstream

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
        '--pipelines-dir',
        help='Override path to the pipelines.',
    )

    return parser


def main(args=None):
    parser = arg_parser()
    args, remaining = parser.parse_known_args(args)
    log.debug(args)

    # Resolve the requested pipeline
    pipelines_home = jetstream.settings['pipelines']['home'].get()
    if not pipelines_home or not os.path.isdir(pipelines_home):
        raise ValueError('Pipelines are not configured. See the "pipelines" '
                         'section of the settings. ')

    pipeline, constants = jetstream.pipelines.get(args.name)

    # Add the pipeline constants to the kvargs
    if constants:
        remaining.extend(['--file:constants', constants])

    # If the pipeline has a bin directory, prepend the env PATH and make
    # the directory available via env variable PIPELINE_BIN
    if pipeline.bin:
        bin_path = pipeline.path.joinpath(pipeline.manifest['bin'])
        os.environ['PIPELINE_BIN'] = str(bin_path)
        os.environ['PATH'] = f'{bin_path}:{os.environ["PATH"]}'

    context = jetstream.context(
        project=args.project,
        kvargs=remaining,
        separator=args.separator
    )

    workflow = jetstream.render_template(
        path=pipeline.main,
        context=context
    )

    try:
        project = jetstream.Project(args.project)
    except jetstream.NotAProject:
        project = None

    if project:
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
