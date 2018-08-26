"""Interact with jetstream projects."""
import sys
import argparse
import logging
import jetstream

log = logging.getLogger(__name__)


def arg_parser(actions=None):
    parser = argparse.ArgumentParser(
        prog='jetstream project',
        description=__doc__
    )

    parser.add_argument('action', choices=actions,
                        help='action name')

    parser.add_argument('args', nargs=argparse.REMAINDER,
                        help='remaining args are passed to the chosen action')

    return parser


def config_arg_parser():
    """Argument parser for the config action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project config',
        description='Summarize project configuration data'
    )

    parser.add_argument('path', nargs='?', default='.',
                        help='Path to a Jetstream project')

    parser.add_argument('--format', choices=['yaml', 'json'], default='yaml')

    parser.add_argument('--json', dest='format',
                        action='store_const', const='json',
                        help='Output JSON')

    parser.add_argument('--yaml', dest='format',
                        action='store_const', const='yaml',
                        help='Output YAML')

    parser.add_argument('--parsable', action='store_true', default=False,
                        help='Parsable output')

    return parser


def tasks_arg_parser():
    """Argument parser for the data action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project tasks',
        description='Summarize project tasks. If `task_id` is not given, all '
                    'task ids will be listed.'
    )

    parser.add_argument('task_id', nargs='*',
                        help='Task ID to summarize')

    parser.add_argument('--path', default='.',
                        help='Jetstream project path')

    return parser


def history_arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream project history',
        description='Records are saved for every workflow that has been run on'
                    'a project. This command lists the run ids in a project.'
    )

    parser.add_argument('path', nargs='?', default='.',
                        help='Path to a Jetstream project')

    return parser


def config(args=None):
    parser =config_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project(args.path)

    if args.format == 'json':
        if args.parsable:
            print(jetstream.utils.json.dumps(p.config))
        else:
            print(jetstream.utils.json.dumps(p.config, indent=4))
    elif args.format == 'yaml':
        jetstream.utils.yaml.dump(p.config, stream=sys.stdout)


def tasks(args=None):
    parser = tasks_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project(args.path)
    wf = p.workflow()
    tasks = {t.tid: t for t in wf.tasks(objs=True)}

    if args.task_id:
        for task_id in args.task_id:
            print(jetstream.utils.yaml_dumps(tasks[task_id].serialize()))
    else:
        print('\t'.join((
            'status',
            'task_id',
            'task_name',
            'start',
            'end'
        )))
        for t in tasks.values():
            print('\t'.join((
                t.state['status'],
                t.tid,
                str(t.directives.get('name')),
                str(t.state['start']),
                str(t.state['end'])
            )))


def history(args=None):
    parser = history_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project(args.path)

    for r in p.history(paths=True):
        r = jetstream.utils.yaml_load(r)
        print(r['id'], r['datetime'])


def main(args=None):
    actions = {
        'config': config,
        'tasks': tasks,
        'history': history,

    }

    parser = arg_parser(actions=list(actions.keys()))
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))
    actions[args.action](args.args)
