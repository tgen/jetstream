"""Manage the tasks in a project or workflow. """
import argparse
import logging
import jetstream

log = logging.getLogger('jetstream.cli')

TASK_DETAILS = """\
- name: {{ task.name }}
  {% if task.directives.cmd %}
  cmd: |
    {{ task.directives.cmd|indent(4) }}
  directives:
  {% for k, v in task.directives.items() if k != 'cmd' %}
    {{ k }}: {{ v|tojson|safe }}
  {% endfor %}
  {% else %}
  directives:
  {% for k, v in task.directives.items() %}
    {{ k }}: {{ v|tojson|safe }}
  {% endfor %}
  {% endif %}
  state:
  {% for k, v in task.state.items() %}
    {{ k }}: {{ v|tojson|safe }}
  {% endfor %}
  {% if logs %}
  logs: |
    {{ logs|indent(4) }}
  {% endif %}
"""

def add_arguments(p):
    p.add_argument(
        'task_names',
        nargs='*',
        default=[],
        help='task names (or patterns) to target'
    )

    filters = p.add_argument_group('task selection')

    filters.add_argument(
        '-s', '--status',
        action='append',
        default=[],
        help='only include tasks with the given status (can be used multiple '
             'times)'
    )

    filters.add_argument(
        '-f', '--format',
        default='glob',
        choices=['exact', 'glob', 'regex'],
        help='change how task names are matched [%(default)s]'
    )

    filters.add_argument(
        '-a', '--ancestors',
        action='store_true',
        help='also select any ancestors of the tasks'
    )

    filters.add_argument(
        '-d', '--descendants',
        action='store_true',
        help='also select any descendants of the tasks'
    )

    actions = p.add_argument_group('actions')

    actions.add_argument(
        '-v', '--verbose',
        action='store_const',
        const='verbose',
        dest='action',
        help='show detailed task info'
    )

    actions.add_argument(
        '--remove',
        action='store_const',
        const='remove',
        dest='action',
        help='remove the targetted tasks from the workflow'
    )

    actions.add_argument(
        '--reset',
        action='store_const',
        const='reset',
        dest='action',
        help='reset the targetted tasks'
    )

    actions.add_argument(
        '--complete',
        action='store_const',
        const='complete',
        dest='action',
        help='complete the targetted tasks'
    )

    actions.add_argument(
        '--fail',
        action='store_const',
        const='fail',
        dest='action',
        help='fail the targetted tasks'
    )

    p.add_argument(
        '-w', '--workflow',
        default=None,
        help='path to a Jetstream workflow file'
    )

    p.add_argument(
        '--no-logs',
        action='store_false',
        dest='include_logs',
        help='dont fetch log files when showing task details (see -v/--verbose)'
    )

    p.add_argument(
        '--action',
        default=None,
        help=argparse.SUPPRESS
    )

    return p


def _iter_targetted_tasks(args, wf):
    """Tasks can be selected by name or pattern, this generator yields any
    tasks covered by the args in the given workflow. This also yields
    ancestors/descendants if those args are set"""
    if args.task_names:
        # We need to search the workflow and possibly query the graph
        if args.descendants or args.ancestors:
            wf.reload_graph()

        for task in _search_workflow(args, wf):
            yield task

            if args.descendants:
                for task in wf.graph.descendants(task):
                    yield task

            if args.ancestors:
                for task in wf.graph.ancestors(task):
                    yield task
    else:
        # Just report on all tasks in the workflow
        for task in wf.tasks.values():
            yield task


def _search_workflow(args, wf):
    """Finds tasks in a workflow by their name, pattern, or status"""
    for task_name in args.task_names:
        if args.format == 'exact':
            yield wf[task_name]
        else:
            for task in wf.find(task_name, style=args.format):
                yield task


def get_details(task, include_logs=True, env=jetstream.templates.environment()):
    if not include_logs or task.is_new():
        logs = None
    else:
        stdout_path = task.state.get('stdout_path')
        if stdout_path is None:
            logs = f'Could not determine log file path.'
        else:
            try:
                with open(stdout_path, 'r') as fp:
                    logs = fp.read()
            except FileNotFoundError:
                logs = f'Could not find log file: {stdout_path}'

    template = env.from_string(TASK_DETAILS)
    final = template.render(task=task, logs=logs)
    return final


def get_summary(task, fields):
    d = task.to_dict()
    values = [jetstream.utils.dict_lookup_dot_notation(d, v) for v in fields]
    return '\t'.join(values)


def main(args):
    log.debug(f'{__name__} {args}')

    if args.workflow:
        log.debug(f'Workflow given by arguments, loading {args.workflow}')
        workflow = jetstream.load_workflow(args.workflow)
    elif args.project:
        log.debug(f'Workflow not given, loading project {args.project}')
        workflow = args.project.load_workflow()
    else:
        err = 'No workflow found. Must be run inside a project, or use ' \
              'options -p/--project -w/--workflow'
        raise FileNotFoundError(err)

    tasks = list(_iter_targetted_tasks(args, workflow))

    if args.status:
        tasks = [t for t in tasks if t.status in args.status]

    if args.action == 'verbose':
        for task in tasks:
            print(get_details(task, include_logs=args.include_logs))
    elif args.action == 'remove':
        for task in tasks:
            workflow.pop(task.name)
        workflow.save()
    elif args.action == 'reset':
        for task in tasks:
            task.reset()
        workflow.save()
    elif args.action == 'complete':
        for task in tasks:
            task.complete()
        workflow.save()
    elif args.action == 'fail':
        for task in tasks:
            task.fail()
        workflow.save()
    else:
        for task in tasks:
            print(task)




