"""Run a template, module, or workflow.
"""
import argparse
import logging
import jetstream
from jetstream.cli import add_config_options_to_parser

log = logging.getLogger('jetstream.cli')


def arg_parser(parser):
    parser.add_argument(
        'path',
        help='path to a template, module, workflow file or pipeline name'
    )

    parser.add_argument(
        '-o', '--out',
        help='path to save the workflow progress [automatically set if working '
             'with a project]'
    )

    parser.add_argument(
        '-b', '--build-only',
        action='store_true',
        help='just render the template, build the workflow, and stop'
    )

    parser.add_argument(
        '-r', '--render-only',
        action='store_true',
        help='just render the template and stop'
    )

    parser.add_argument(
        '--backend',
        choices=jetstream.settings['backends'].get(dict),
        default=jetstream.settings['backend'].get(str),
        help='runner backend name used for executing tasks [%(default)s]'
    )

    parser.add_argument(
        '--format',
        choices=['template', 'module', 'workflow'],
        help='workflow format [inferred from the extension of the path]'
    )

    parser.add_argument(
        '--reset-method',
        choices=['retry', 'resume', 'reset'],
        default='retry',
        help='controls which tasks are reset prior to starting the run - '
             '"retry": pending and failed, "resume": pending, or "reset": '
             'all [%(default)s]'
    )

    parser.add_argument(
        '--existing-workflow',
        help='path to an existing workflow file that will be merged into run '
             '[automatically set if working with a project]'
    )

    parser.add_argument(
        '--template-dir',
        action='append',
        default=[],
        help='directory to add to search path for loading templates, this can '
             'be used multiple times [automatically set with pipelines command]'
    )

    parser.add_argument(
        '--pipeline',
        help=argparse.SUPPRESS
    )

    parser.add_argument(
        '--runner',
        help=argparse.SUPPRESS
    )

    add_config_options_to_parser(parser)
    return parser


def _resolve_format(args):
    """Set the workflow format arg based on the path extension"""

    # TODO This is not consistent with the file extension conventions
    # set in the workflow.save() method. Workflows and templates are both
    # allowed flexible file extensions. We could sniff the file contents, but
    # that creates issues with pipes.
    if args.format is None:
        if args.path.endswith('.pickle'):
            return 'workflow'
        elif args.path.endswith('.py'):
            return 'module'
        else:
            return 'template'
    else:
        return args.format


def render_only(args):
    render = jetstream.templates.render_template(
        path=args.path,
        project=args.project,
        pipeline=args.pipeline,
        command_args=args.config,
        search_path=args.template_dir
    )
    if args.out:
        with open(args.out, 'w') as fp:
            fp.write(render)
    else:
        print(render)


def build_only(args):
    wf = jetstream.templates.build_template(
        path=args.path,
        project=args.project,
        pipeline=args.pipeline,
        command_args=args.config,
        search_path=args.template_dir
    )
    wf.graph()
    if args.out:
        wf.save(args.out)
    else:
        print(wf)


def from_module(args):
    raise NotImplementedError


def from_template(args):
    return jetstream.templates.build_template(
        path=args.path,
        project=args.project,
        pipeline=args.pipeline,
        command_args=args.config
    )


def from_workflow(args):
    return jetstream.load_workflow(args.path)


def check_for_failures(tasks):
    fails = []
    skips = []
    for task in tasks:
        if task.is_failed():
            if task.is_skipped():
                skips.append(task)
            else:
                fails.append(task)

    if fails or skips:
        log.info(f'{len(fails)+len(skips)} total tasks failed:')
        for i in range(5):
            try:
                task = fails.pop(-1)
                log.info(f'{task} {task.state.get("label")}')
            except IndexError:
                break
        else:
            if fails:
                log.info(f'...and {len(fails)} more failed!')

        if skips:
            log.info(f'{len(skips)} tasks were skipped due to '
                     f'failed dependencies')

        raise RuntimeError('Some tasks failed')


def run(args):
    """Given some new workflow and the command args, this will:
    1) load any existing workflow either from specific args or project
    2) mash the existing workflow with the new workflow
    3) run the combined workflow
    4) raise a RuntimeError if any tasks failed
    """
    args.runner = jetstream.Runner(backend=args.backend)

    format = _resolve_format(args)

    if args.render_only:
        return render_only(args)
    elif args.build_only:
        return build_only(args)
    elif format == 'template':
        wf = from_template(args)
    elif format == 'module':
        wf = from_module(args)
    elif format == 'workflow':
        wf = from_workflow(args)
    else:
        msg = f'Invalid argument configuration: {args}'
        raise ValueError(msg)

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

    try:
        args.runner.start(
            workflow=wf,
            project=args.project
        )
    finally:
        check_for_failures(wf.tasks.values())



def main(args):
    log.debug(f'{__name__} {args}')

    try:
        if args.project:
            with args.project.lock:
                run(args)
        else:
            run(args)
    except RuntimeError:
        log.error('There were errors during the run.')
        raise SystemExit(1)
    except TimeoutError:
        log.error('Failed to acquire project lock, there may be a run pending.')
        raise SystemExit(1)

