"""Draw the network graph for a workflow.

This command will load a finalized workflow and create an image of the network
graph

This requires several dependencies that are not automatically installed with
Jetstream. TODO: Add optional install params to setup.py?

"""
import logging
import jetstream

log = logging.getLogger(__name__)


def arg_parser(p):
    p.add_argument(
        'out',
        help='Path to save the image file'
    )

    p.add_argument(
        '-w', '--workflow',
        help='Path to a Jetstream workflow file'
    )


def main(args):
    log.debug(f'{__name__} {args}')

    if args.workflow:
        args.workflow = jetstream.load_workflow(args.workflow)
    else:
        if args.project is None:
            raise ValueError('No workflow given and not working in a project')
        args.workflow = args.project.workflow

    workflow = jetstream.load_workflow(args.workflow)
    workflow.draw(filename=args.out)
