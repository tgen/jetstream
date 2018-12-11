"""Interact with jetstream projects. View tasks, run history, or project data.
This option requires a subcommand.
"""
import sys
import argparse
import logging
import jetstream
import textwrap
from jetstream.cli import shared

log = logging.getLogger(__name__)


def arg_parser():
    parent = argparse.ArgumentParser(add_help=False)

    parent.add_argument(
        '--project',
        default=None,
        help='If the cwd is a project, it will be loaded automatically. '
             'Otherwise, a path to a project can be specified.'
    )

    parser = argparse.ArgumentParser(
        prog='jetstream project',
        add_help=False,
        description=__doc__,
        parents=[parent,]
    )

    parser.add_argument(
        '-h',
        '--help',
        dest='subcommand',
        action='store_const',
        const='help'
    )

    subparsers = parser.add_subparsers(
        title='subcommands',
        dest='subcommand',
        help=argparse.SUPPRESS
    )

    subparsers.add_parser(
        name='help',
        description='Show detailed help for the project command',
    )

    variables = subparsers.add_parser(
        name='variables',
        description='Show project config variables',
        parents=[parent,]
    )

    variables.add_argument(
        '--format',
        choices=['yaml', 'json'],
        default='yaml')

    variables.add_argument(
        '--json',
        dest='format',
        action='store_const',
        const='json',
        help='Output JSON')

    variables.add_argument(
        '--yaml',
        dest='format',
        action='store_const',
        const='yaml',
        help='Output YAML'
    )

    variables.add_argument(
        '--minified',
        action='store_true',
        default=False,
        help='Minified output instead of human-friendly'
    )

    tasks = subparsers.add_parser(
        name='tasks',
        description='Summarize project tasks. If `task_id` is not given, all '
                    'tasks will be listed.',
        parents = [parent, ]
    )

    tasks.add_argument(
        'task_name',
        nargs='*',
        help='Task name(s) to summarize'
    )

    task_remove = subparsers.add_parser(
        name='remove_tasks',
        description='Warning - Experimental feature! '
                    'Remove tasks from the current project workflow',
        parents = [parent, ]
    )

    task_remove.add_argument(
        'task_name',
        nargs='+',
        help='Task name(s) to remove. This will search for tasks with a '
             'name matching the pattern and remove them from the workflow.'
    )

    task_remove.add_argument(
        '-d', '--descendants',
        action='store_true',
        default=False,
        help='Also remove any descendants of the task'
    )

    task_remove.add_argument(
        '-f', '--force',
        action='store_true',
        default=False,
        help='Ignore dependency errors. Warning - This may corrupt the '
             'workflow causing it to be unloadable.'
    )

    task_reset = subparsers.add_parser(
        name='reset_tasks',
        description='Warning - Experimental feature! '
                    'Reset tasks in the current project workflow',
        parents=[parent, ]
    )

    task_reset.add_argument(
        'task_name',
        nargs='+',
        help='Task name(s) to reset'
    )

    task_complete = subparsers.add_parser(
        name='complete_tasks',
        description='Warning - Experimental feature! '
                    'Complete tasks in the current project workflow',
        parents=[parent, ]
    )

    task_complete.add_argument(
        'task_name',
        nargs='+',
        help='Task name(s) to complete'
    )

    task_fail = subparsers.add_parser(
        name='fail_tasks',
        description='Warning - Experimental feature! '
                    'Fail tasks in the current project workflow',
        parents=[parent, ]
    )

    task_fail.add_argument(
        'task_name',
        nargs='+',
        help='Task name(s) to fail'
    )

    history = subparsers.add_parser(
        name='history',
        description='Records are saved for every run that has been started on'
                    'a project. This command lists the run ids in a project.',
        parents=[parent, ]
    )

    history.add_argument(
        'run_id',
        nargs='*',
        help='Run ID to summarize'
    )

    return parser


class Subcommands:
    @staticmethod
    def variables(args):
        p = jetstream.Project(path=args.project)

        if args.format == 'json':
            if args.minified:
                print(jetstream.utils.json_dumps(p.config))
            else:
                print(jetstream.utils.json_dumps(p.config, indent=4))
        elif args.format == 'yaml':
            jetstream.utils.yaml_dump(p.config, stream=sys.stdout)

    @staticmethod
    def tasks(args=None):
        p = jetstream.Project(path=args.project)
        wf = p.workflow()

        if args.task_name:
            for task_name in args.task_name:
                tids = wf.find(task_name)
                for tid in tids:
                    task = wf.get_task(tid)
                    print(jetstream.utils.yaml_dumps(task.serialize()))
        else:
            tasks = {t.tid: t for t in wf.tasks(objs=True)}

            print('\t'.join(('task_id', 'identity', 'status',)))
            for t in tasks.values():
                print('\t'.join((t.tid, t.identity, t.status,)))

    @staticmethod
    def remove_tasks(args=None):
        p = jetstream.Project(path=args.project)
        wf = p.workflow()

        for name in args.task_name:
            wf.remove_task(name, force=args.force, descendants=args.descendants)

        wf.save()

    @staticmethod
    def reset_tasks(args=None):
        p = jetstream.Project(path=args.project)
        wf = p.workflow()

        for name in args.task_name:
            task_ids = wf.find(name)
            for tid in task_ids:
                wf.get_task(tid).reset()

        wf.save()

    @staticmethod
    def complete_tasks(args=None):
        p = jetstream.Project(path=args.project)
        wf = p.workflow()

        for name in args.task_name:
            task_ids = wf.find(name)
            for tid in task_ids:
                wf.get_task(tid).complete()

        wf.save()

    @staticmethod
    def fail_tasks(args=None):
        p = jetstream.Project(path=args.project)
        wf = p.workflow()

        for name in args.task_name:
            task_ids = wf.find(name)
            for tid in task_ids:
                wf.get_task(tid).fail()

        wf.save()

    @staticmethod
    def history(args=None):
        p = jetstream.Project(path=args.project)

        for r in p.history(paths=True):
            r = jetstream.utils.load_yaml(r)
            print(r['id'], r['datetime'])

    @staticmethod
    def help(parser):
        subparsers = [action for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)]
        print(parser.format_help())
        print('Available subcommands:\n')
        for s in subparsers:
            for choice, subparser in s.choices.items():
                if choice == 'help':
                    continue
                else:
                    desc = textwrap.wrap(subparser.description, width=70)
                    desc = '\n'.join(desc)
                    text = f'{choice}:\n\n{textwrap.indent(desc, "  ")}\n'
                    print(textwrap.indent(text, "  "))


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    if args.project is None:
        try:
            jetstream.Project()
        except jetstream.NotAProject:
            raise ValueError(
                'This command requires a jetstream project. Run '
                'inside a jetstream project or use -p/--project'
            )

    if args.subcommand == 'help':
        Subcommands.help(parser)
    else:
        getattr(Subcommands, args.subcommand)(args)

