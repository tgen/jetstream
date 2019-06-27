"""Manage the tasks in a workflow"""
import argparse
import logging
import jetstream

log = logging.getLogger('jetstream.cli')


def arg_parser(p):
    p.add_argument(
        '-w', '--workflow',
        default=None,
        help='path to a Jetstream workflow file'
    )

    p.add_argument(
        '-n', '--name',
        nargs='+',
        default=[],
        help='task name(s)'
    )

    p.add_argument(
        '-r', '--regex',
        nargs='+',
        default=[],
        help='task name pattern(s) - match task names with regex'
    )

    shared = argparse.ArgumentParser(add_help=False)

    shared.add_argument(
        '-w', '--workflow',
        default=None,
        help='path to a Jetstream workflow file'
    )

    shared.add_argument(
        '-n', '--name',
        nargs='+',
        default=[],
        help='task name(s)'
    )

    shared.add_argument(
        '-r', '--regex',
        nargs='+',
        default=[],
        help='task name pattern(s) - match task names with regex'
    )

    subparsers = p.add_subparsers(
        dest='subcommand'
    )

    # Task detailed view
    detailed_tasks = subparsers.add_parser(
        name='details',
        parents=[shared,],
        description='See a more detailed view of tasks'
    )

    detailed_tasks.set_defaults(func=details)

    # Removing tasks
    remove_tasks = subparsers.add_parser(
        name='remove',
        parents=[shared, ],
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
        name='reset',
        parents=[shared, ],
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
        name='complete',
        parents=[shared, ],
        description='Complete tasks in the workflow'
    )

    complete_tasks.set_defaults(func=complete)

    # Fail tasks
    fail_tasks = subparsers.add_parser(
        name='fail',
        parents=[shared, ],
        description='Fail tasks in the workflow'
    )

    fail_tasks.set_defaults(func=fail)


def _iter_selected_tasks(args, wf):
    """Tasks can be selected by name or pattern, this generator yields any
    tasks covered by the args in the given workflow. """
    for task_name in args.name:
        try:
            yield wf.tasks[task_name]
        except KeyError as e:
            log.error(f'"{task_name}" not found in workflow')


    for pat in args.regex:
        for t in wf.find(pat):
            yield t


def details(args):
    """Long format details about the tasks"""
    log.debug(f'{__name__} {args}')
    wf = load_workflow(args)

    if args.name or args.regex:
        tasks = set()
        for task in _iter_selected_tasks(args, wf):
            tasks.add(task)
    else:
        tasks = wf.tasks.values()

    for t in tasks:
        print(t.name)
        print(f'Directives:\n{jetstream.utils.yaml_dumps(t.directives)}')
        print(f'State:\n{jetstream.utils.yaml_dumps(t.state)}')
        print('Logs:')
        stdout_path = t.state.get('stdout_path') or t.directives.get('stdout')
        if stdout_path is None:
            print(f'Could not determine the log path')
        else:
            try:
                with open(stdout_path, 'r') as fp:
                    print(fp.read())
            except FileNotFoundError:
                print(f'Could not find log file: {stdout_path}')


def remove(args):
    log.debug(f'{__name__} {args}')
    wf = load_workflow(args)

    msg = '{} has successors that would be orphaned, remove them first' \
          'or use --force option.'

    if args.force:
        for task in _iter_selected_tasks(args, wf):
            log.info(f'Removing: {task}')
            wf.pop(task.name)
    else:
        graph = wf.graph()

        for task in _iter_selected_tasks(args, wf):
            has_deps = next(iter(graph.successors(task)))
            if has_deps:
                raise ValueError(msg.format(task))
            else:
                log.info(f'Removing: {task}')
                wf.pop(task.name)

    wf.save()


def reset(args):
    log.debug(f'{__name__} {args}')
    wf = load_workflow(args)

    for task in _iter_selected_tasks(args, wf):
        log.info(f'Resetting: {task}')
        task.reset()

    wf.save()


def complete(args):
    log.debug(f'{__name__} {args}')
    wf = load_workflow(args)

    for task in _iter_selected_tasks(args, wf):
        log.info(f'Completing: {task}')
        task.complete()

    wf.save()


def fail(args):
    log.debug(f'{__name__} {args}')
    wf = load_workflow(args)

    for task in _iter_selected_tasks(args, wf):
        log.info(f'Failing: {task}')
        task.fail(force=True)

    wf.save()


def summary(args):
    log.debug(f'{__name__} {args}')
    wf = load_workflow(args)

    if args.name or args.regex:
        tasks = set()
        for task in _iter_selected_tasks(args, wf):
            tasks.add(task)
    else:
        tasks = list(wf)

    f = jetstream.settings['tasks']['summary_fields'].get(list)
    print('\t'.join(f))

    for t in tasks:
        d = t.to_dict()
        values = [jetstream.utils.dict_lookup_dot_notation(d, v) for v in f]
        print('\t'.join(values))


def load_workflow(args):
    if args.workflow:
        log.debug(f'Workflow given by arguments, loading {args.workflow}')
        workflow = jetstream.load_workflow(args.workflow)
    elif args.project:
        log.debug(f'Workflow not given by arguments, loading project {args.project}')
        workflow = args.project.load_workflow()
    else:
        err = 'No workflow given. Must be run inside a project, or use ' \
              '--project --workflow arguments'
        raise FileNotFoundError(err)
    return workflow


def main(args):
    log.debug(f'{__name__} {args}')
    summary(args)
