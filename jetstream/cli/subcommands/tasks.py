"""Manage the tasks in a workflow"""
import logging
import jetstream

log = logging.getLogger('jetstream.cli')


def arg_parser(p):
    p.add_argument(
        '-w', '--workflow',
        help='path to a Jetstream workflow file'
    )

    p.add_argument(
        '-n', '--name',
        nargs='+',
        help='task name(s)'
    )

    p.add_argument(
        '-r', '--regex',
        nargs='+',
        help='task name pattern(s) - match task names with regex'
    )

    subparsers = p.add_subparsers(
        dest='subcommand'
    )

    # Task detailed view
    detailed_tasks = subparsers.add_parser(
        name='details',
        description='See a more detailed view of tasks'
    )

    detailed_tasks.set_defaults(func=details)

    # Removing tasks
    remove_tasks = subparsers.add_parser(
        name='remove',
        description='Remove tasks from the current project workflow'
    )

    remove_tasks.set_defaults(func=remove)


    remove_tasks.add_argument(
        '-d', '--descendants',
        action='store_true',
        help='also remove any descendants of the task'
    )

    remove_tasks.add_argument(
        '-f', '--force',
        action='store_true',
        help='Ignore dependency errors. Caution! This may leave tasks with '
             'missing dependencies and result in a corrupted workflow.'
    )

    # Reset tasks
    reset_tasks = subparsers.add_parser(
        name='reset_tasks',
        aliases=['reset', ],
        description='Reset tasks in the workflow. '
    )

    reset_tasks.set_defaults(func=reset)

    reset_tasks.add_argument(
        '-a', '--ancestors',
        action='store_true',
        help='also reset any ancestors of the task'
    )

    # Complete tasks
    complete_tasks = subparsers.add_parser(
        name='complete_tasks',
        aliases=['complete', ],
        description='Complete tasks in the workflow'
    )

    complete_tasks.set_defaults(func=complete)

    # Fail tasks
    fail_tasks = subparsers.add_parser(
        name='fail_tasks',
        aliases=['fail', ],
        description='Fail tasks in the workflow'
    )

    fail_tasks.set_defaults(func=fail)


def summary(args):
    wf = args.workflow

    tasks = set()
    if args.name or args.regex:
        for task_name in args.name:
            task = wf.get_task(task_name)
            tasks.add(task)
        for pat in args.pattern:
            for t in wf.find(pat, objs=True):
                tasks.add(t)
    else:
        tasks = list(wf)

    f = jetstream.settings['tasks']['summary_fields'].get(list)
    print('\t'.join(f))

    for t in tasks:
        d = t.to_dict()
        values = [jetstream.utils.dict_lookup_dot_notation(d, v) for v in f]
        print('\t'.join(values))


def details(args):
    wf = args.workflow

    tasks = set()
    if args.name or args.regex:
        for task_name in args.name:
            task = wf.get_task(task_name)
            tasks.add(task)
        for pat in args.regex:
            for t in wf.find(pat, objs=True):
                tasks.add(t)
    else:
        tasks = wf.list_tasks()

    for t in tasks:
        # TODO Finish this
        print(t.to_dict())


def remove(args):
    wf = args.workflow
    msg = '{} has successors that would be orphaned, remove them first' \
          'or use --force option.'

    if args.force:
        for name in args.task_name:
            wf.pop(name)

        for pat in args.regex:
            for task in wf.find(pat, fallback=[]):
                wf.pop(task.name)
    else:
        graph = wf.graph()

        for name in args.task_name:
            task = wf[name]
            has_deps = next(iter(graph.successors(task)))
            if has_deps:
                raise ValueError(msg.format(task))
            else:
                wf.pop(name)

        for pat in args.regex:
            for task in wf.find(pat, fallback=[]):
                has_deps = next(iter(graph.successors(task)))
                if has_deps:
                    raise ValueError(msg.format(task))
                else:
                    wf.pop(task.name)
    wf.save()


def reset(args):
    wf = args.workflow
    for name in args.task_name:
        task_ids = wf.find(name)
        for name in task_ids:
            wf[name].reset()
    wf.save()


def complete(args):
    wf = args.workflow
    for name in args.task_name:
        task_ids = wf.find(name)
        for name in task_ids:
            wf[name].complete()
    wf.save()


def fail(args):
    wf = args.workflow
    for name in args.task_name:
        task_ids = wf.find(name)
        for name in task_ids:
            wf[name].fail()
    wf.save()


def main(args):
    log.debug(f'{__name__} {args}')

    if args.workflow:
        args.workflow = jetstream.load_workflow(args.workflow)
    elif args.project:
        args.workflow = args.project.load_workflow()
    else:
        err = 'No workflow given. Must be run inside a project, or use ' \
              '--project --workflow arguments'
        raise FileNotFoundError(err)


    if args.func is main:
        summary(args)

