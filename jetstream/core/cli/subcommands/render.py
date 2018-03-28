"""Command line utility for rendering workflow templates with project data.
Multiple data files can be used during the render process in a "last value
wins" pattern. Jinja2 will raise errors when trying to access nested values
for which a parent value is undefined. However, undefined leaves in the data
will only generate a warning unless "--strict" is used.
"""
import argparse
import json
import yaml
import logging
from jetstream.core.workflows.build import render_template, parse_nodes

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('template', help='Path to a workflow template')

    parser.add_argument('data', nargs='*',
                        help='Path to a run data file(s). Used to render'
                             'a workflow template.')

    parser.add_argument('--strict', action='store_true', default=False)

    parser.add_argument('--render-only', action='store_true', default=False)

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    data = {}
    if args.data:
        for path in args.data:
            with open(path, 'r') as fp:
                raw = fp.read()

            if path.endswith('.yaml'):
                data.update(yaml.load(raw))
            elif path.endswith('.json'):
                data.update(json.loads(raw))
            else:
                # TODO allow explicit override of file types
                raise ValueError('Unrecognized run data format {}'.format(path))

    with open(args.template, 'r') as fp:
        template = fp.read()
    log.debug('Template:\n{}'.format(template))

    render = render_template(template, data, strict=args.strict)
    log.debug('Render:\n{}'.format(render))

    if args.render_only:
        print(render)
    else:
        parsed = yaml.load(render)
        log.debug('Parsed:\n{}'.format(parsed))

        print(parse_nodes(parsed))

