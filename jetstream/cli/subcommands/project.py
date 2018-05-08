""" This module contains the cli interface code for the project data utility."""
import os
import sys
import argparse
import logging
import jetstream

log = logging.getLogger(__name__)


def arg_parser(actions):
    parser = argparse.ArgumentParser(
        prog='jetstream project',
        description=__doc__
    )

    parser.add_argument('action', choices=actions)

    parser.add_argument('args', nargs=argparse.REMAINDER)

    return parser


def init_arg_parser():
    """arg parser for the init action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project init',
    )

    parser.add_argument('path', nargs='?', default='.')

    return parser


def config_arg_parser():
    """arg parser for the config action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project config',
    )

    parser.add_argument('path', nargs='?', default='.')

    parser.add_argument('--format', choices=['yaml', 'json'], default='json')

    parser.add_argument('--json', dest='format',
                        action='store_const', const='json')

    parser.add_argument('--yaml', dest='format',
                        action='store_const', const='yaml')

    parser.add_argument('--pretty', action='store_true', default=False)

    return parser


def task_summary_arg_parser():
    """arg parser for the data action"""
    parser = argparse.ArgumentParser(
        prog='jetstream project task_summary',
    )

    parser.add_argument('task_id', nargs='*')

    parser.add_argument('--path', default='.')

    parser.add_argument('-r', '--run', default='latest')

    return parser


def runs_arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream project runs'
    )

    parser.add_argument('path', nargs='?', default='.')

    return parser


def init(args=None):
    parser = init_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    jetstream.project.init(args.path)


def config(args=None):
    parser =config_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project(args.path)

    if args.format == 'json':
        if args.pretty:
            print(jetstream.utils.json.dumps(p.config, indent=4))
        else:
            print(jetstream.utils.json.dumps(p.config))
    elif args.format == 'yaml':
        jetstream.utils.yaml.dump(p.config, stream=sys.stdout)


def samples(args=None):
    parser = config_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project(args.path)

    if args.format == 'json':
        if args.pretty:
            print(jetstream.utils.json.dumps(p.samples(), indent=4))
        else:
            print(jetstream.utils.json.dumps(p.samples()))
    elif args.format == 'yaml':
        jetstream.utils.yaml.dump(p.samples(), stream=sys.stdout)


def task_summary(args=None):
    parser = task_summary_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project(args.path)

    if args.run == 'latest':
        run_id = p.latest_run()
    else:
        if not args.run in p.runs():
            raise ValueError("Run ID {} not found in project.".format(args.run))
        else:
            run_id = args.run

    run_path = os.path.join(p.path, jetstream.project_index, run_id)
    workflow_path = os.path.join(run_path, 'workflow.yaml')
    wf = jetstream.workflows.load(workflow_path)
    tasks = dict(wf.nodes(data=True))

    if args.task_id:
        for task_id in args.task_id:
            task = tasks[task_id]
            print(jetstream.utils.task_summary(task_id, task))
    else:
        for t in tasks.keys():
            print(t)


def runs(args=None):
    parser = runs_arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project(args.path)

    for r in p.runs():
        try:
            created = os.path.join(p.path, p.index_path, r, 'created.yaml')
            with open(created, 'r') as fp:
                print(fp.read())
        except Exception as e:
            log.exception(e)



def main(args=None):
    actions = {
        'init': init,
        'config': config,
        'samples': samples,
        'task_summary': task_summary,
        'runs': runs,

    }

    parser = arg_parser(actions=list(actions.keys()))
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))
    actions[args.action](args.args)
