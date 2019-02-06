"""Render a template and print the results"""
import logging
from jetstream.templates import context, render_template

log = logging.getLogger(__name__)


def arg_parser(parser):
    parser.add_argument(
        'path',
        help='Path to a workflow template'
    )



def main(args):
    log.debug(f'{__name__} {args}')

    c = context(
        project=args.project,
        command_args=args.kvargs
    )

    render = render_template(args.path, c, render_only=True)
    print(render)
