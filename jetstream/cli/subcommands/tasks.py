"""Manage the tasks in a workflow"""
import logging
import sys
import jetstream

log = logging.getLogger(__name__)


def arg_parser(p):
    p.add_argument(
        '-t', '--tasks',
        nargs='+',
        help='Task name(s) to summarize'
    )

    p.add_argument(
        '-w', '--workflow',
        help='Path to a Jetstream workflow file'
    )

    subparsers = p.add_subparsers(
        dest='subcommand'
    )

    remove_tasks = subparsers.add_parser(
        name='remove_tasks',
        aliases=['remove',],
        description='Remove tasks from the current project workflow'
                    'Warning - Experimental feature! If you '
                    'find this feature useful, have problems or suggestions, '
                    'please submit an issue on Github.'
    )

    remove_tasks.set_defaults(func=remove)

    remove_tasks.add_argument(
        'task_name',
        nargs='+',
        help='Task name(s) to remove. This will search for tasks with a '
             'name matching the pattern and remove them from the workflow.'
    )

    remove_tasks.add_argument(
        '-d', '--descendants',
        action='store_true',
        default=False,
        help='Also remove any descendants of the task'
    )

    remove_tasks.add_argument(
        '-f', '--force',
        action='store_true',
        default=False,
        help='Ignore dependency errors. Warning - This may corrupt the '
             'workflow causing it to be unloadable.'
    )

    reset_tasks = subparsers.add_parser(
        name='reset_tasks',
        aliases=['reset', ],
        description='Reset tasks in the current project workflow. '
                    'Warning - Experimental feature! If you '
                    'find this feature useful, have problems or suggestions, '
                    'please submit an issue on Github.'
    )

    reset_tasks.set_defaults(func=reset)

    complete_tasks = subparsers.add_parser(
        name='complete_tasks',
        aliases=['complete', ],
        description='Complete tasks in the current project workflow'
                    'Warning - Experimental feature! If you '
                    'find this feature useful, have problems or suggestions, '
                    'please submit an issue on Github.'
    )

    complete_tasks.set_defaults(func=complete)

    fail_tasks = subparsers.add_parser(
        name='fail_tasks',
        aliases=['fail', ],
        description='Fail tasks in the current project workflow'
                    'Warning - Experimental feature! If you '
                    'find this feature useful, have problems or suggestions, '
                    'please submit an issue on Github.'
    )

    fail_tasks.set_defaults(func=fail)


def tasks(args):
    wf = args.workflow

    if args.tasks:
        for task_name in args.tasks:
            tasks = wf.find(task_name, objs=True)
            info = [task.serialize() for task in tasks]
            jetstream.utils.yaml_dump(info, sys.stdout)
    else:
        tasks = {t.tid: t for t in wf.tasks(objs=True)}

        print('\t'.join(('task_id', 'status',)))
        for t in tasks.values():
            print('\t'.join((t.tid, t.status,)))


def remove(args):
    wf = args.workflow

    for name in args.task_name:
        wf.remove_task(name, force=args.force, descendants=args.descendants)

    wf.save()


def reset(args):
    wf = args.workflow

    for name in args.task_name:
        task_ids = wf.find(name)
        for tid in task_ids:
            wf.get_task(tid).reset()

    wf.save()


def complete(args):
    wf = args.workflow

    for name in args.task_name:
        task_ids = wf.find(name)
        for tid in task_ids:
            wf.get_task(tid).complete()

    wf.save()


def fail(args):
    wf = args.workflow

    for name in args.task_name:
        task_ids = wf.find(name)
        for tid in task_ids:
            wf.get_task(tid).fail()

    wf.save()


def main(args):
    log.debug(f'{__name__} {args}')

    if args.workflow:
        args.workflow = jetstream.load_workflow(args.workflow)
    else:
        if args.project is None:
            raise ValueError('No workflow given and not working in a project')
        args.workflow = args.project.workflow

    if args.func is main:
        tasks(args)

