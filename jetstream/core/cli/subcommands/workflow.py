"""Command line utility for one-step render and run of workflow templates"""
import os
import argparse
import logging
from jetstream.core.project import Project
from jetstream.core.workflows.builder import render_template, build_workflow
from jetstream.core.run import run_workflow
from jetstream import utils


log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('template', help='Path to a workflow template')

    parser.add_argument('data', nargs='*',
                        help='Path to a run data file(s). Used to render'
                             'a workflow template.')

    parser.add_argument('--strict', action='store_true', default=False)

    return parser


def render(template, data, strict=True):
    all_data = {}

    for path in data:
        name = os.path.splitext(os.path.basename(path))[0]


        with open(path, 'r') as fp:
            raw = fp.read()

        if path.endswith('.yaml'):
            obj = {name: utils.yaml_loads(raw)}
            all_data.update(obj)
        elif path.endswith('.json'):
            obj = {name: utils.json_loads(raw)}
            all_data.update(obj)
        elif path.endswith(('.csv', '.tsv')):
            obj = {name: utils.table_to_records(path)}
            all_data.update(obj)
        else:
            # TODO allow explicit override of file types
            raise ValueError('Unrecognized run data format {}'.format(path))

    with open(template, 'r') as fp:
        template = fp.read()

    log.debug('Template:\n{}'.format(template))

    r = render_template(template, all_data, strict=strict)
    log.debug('Render:\n{}'.format(r))

    return r


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    # Load the project, ensure we're working in a valid project
    p = Project()

    # Converts a template into a workflow
    rendered_template = render(
        template=args.template,
        strict=args.strict,
        data=args.data
    )

    # Rendered template is a yaml format array of nodes
    # we load this in with the yaml library, then build a
    # workflow from the nodes
    nodes = utils.yaml.load(rendered_template)
    wf = build_workflow(nodes)

    # Now we run the workflow in the project
    run_workflow(wf, p)
