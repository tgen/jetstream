"""Run a template, module, workflow, or pipeline"""
import argparse
import logging
import subprocess
import jetstream
from jetstream.cli.subcommands import run_common_options

log = logging.getLogger('jetstream.cli')


def add_arguments(parser):
    parser.add_argument(
        'file',
        nargs='?',
        help='path to a template, module, workflow, or pipeline'
    )

    parser.add_argument(
        '--format',
        choices=('template', 'module', 'workflow'),
        default='template',
        help='file format [%(default)s]'
    )

    run_common_options.add_arguments(parser)


def task_log_snippets(task, project):
    stdin, stdout, stderr = jetstream.tasks.get_fd_paths(task, project)

    if stdout:
        logs = jetstream.utils.get_snippet(stdout)
        log.info(f'Stdout snippet:\n{logs}')

    if stderr != stdout:
        logs = jetstream.utils.get_snippet(stderr)
        log.info(f'Stderr snippet:\n{logs}')


def check_for_failures(wf, project=None):
    limit = jetstream.settings['runner']['report_failures'].as_number()
    logs = jetstream.settings['runner']['report_failures_logs'].as_number()
    fails = []
    skips = []
    for task in wf.tasks.values():
        if task.is_failed():
            if task.is_skipped():
                skips.append(task)
            else:
                fails.append(task)

    if fails or skips:
        log.info(f'{len(fails) + len(skips)} total tasks failed:')
        remaining = limit
        while remaining:
            try:
                task = fails.pop(-1)
            except IndexError:
                break

            log.info(f'Task failed: "{task.name}": {task.state.get("label")}')
            if logs:
                task_log_snippets(task, project)

            remaining -= 1

        if fails:
            log.info(f'...and {len(fails)} more failed!')

        if skips:
            log.info(f'{len(skips)} tasks were skipped due to '
                     f'failed dependencies')

        return True
    else:
        return False


def template(template, args):
    render = jetstream.templates.render_template(
        template,
        project=args.project,
        pipeline=args.pipeline,
        command_args=args.config,
    )

    if args.render_only:
        if args.out:
            with open(args.out, 'w') as fp:
                print(render, file=fp)
        else:
            print(render)
            return

    wf = jetstream.templates.load_workflow(render)

    if args.build_only:
        if args.out:
            wf.save(args.out)
        log.info('Workflow built successfully!') 
        return

    run(wf, args)


def workflow(args):
    wf = jetstream.load_workflow(args.file)
    run(wf, args)


def module(args):
    raise NotImplementedError


def run(wf, args):
    """Given some new workflow and the command args, this will:
    1) load any existing workflow either from specific args or project
    2) mash the existing workflow with the new workflow
    3) run the combined workflow
    4) raise a RuntimeError if any tasks failed
    """
    args.runner = jetstream.Runner(backend=args.backend)

    if args.existing_workflow:
        ewf = jetstream.load_workflow(args.existing_workflow)
    elif args.project:
        ewf = args.project.load_workflow()
    else:
        ewf = None

    if ewf is not None:
        wf = jetstream.workflows.mash(ewf, wf)
        wf.path = ewf.path

    if args.out:
        wf.path = args.out

    wf.reset(args.reset_method)

    if args.mash_only:
        if args.out:
            wf.save(args.out)
        else:
            wf.save()
        log.info('Mashed workflow built successfully!')
        return

    try:
        args.runner.start(
            workflow=wf,
            project=args.project,
            pipeline=args.pipeline
        )
    finally:
        errs = check_for_failures(wf, args.project)
        if errs:
            log.critical('There were errors during the run.')
            raise SystemExit(1)


def main(args):
    log.debug(f'{__name__} {args}')

    if not (args.pipeline or args.file):
        raise ValueError('Must give file path if pipeline is not set')

    if args.pipeline:
        t = args.pipeline.load_template()
        template(t, args)
    elif args.format == 'template':
        t = jetstream.templates.load_template(path=args.file, *args.search_path)
        template(t, args)
    elif args.format == 'workflow':
        workflow(args)
    elif args.format == 'module':
        module(args)
    else:
        raise ValueError(f'unknown file format {args.format}')
