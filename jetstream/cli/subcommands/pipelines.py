"""Run a pipeline

Pipelines are Jetstream templates that have been documented with version
information and added to the jetstream pipelines directory. This command
allows pipelines to be referenced by name and automatically includes the
pipeline scripts and constants in the run."""
import os
import logging
import jetstream

log = logging.getLogger(__name__)


def arg_parser(parser):
    parser.add_argument(
        'name',
        nargs='?',
        help='Pipeline name'
    )

    parser.add_argument(
        '--pipelines-dir',
        help='Override path to the pipelines.',
    )

    return parser


def main(args=None):
    log.debug(f'{__name__} {args}')

    raise NotImplementedError("TODO: FINISH THIS COMMAND")

    if args.name is None:
        print('WHAT', jetstream.pipelines.ls())
        return

    # Resolve the requested pipeline
    pipelines_home = jetstream.settings['pipelines']['home'].get()
    if not pipelines_home or not os.path.isdir(pipelines_home):
        err = 'Pipelines are not configured. Check the application settings.'
        raise ValueError(err)

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
    except jetstream.ProjectInvalid:
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
