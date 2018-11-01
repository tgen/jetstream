"""Draw the network graph for a workflow.

This command will load a finalized workflow and create an image of the network
graph

This requires several dependencies that are not automatically installed with
Jetstream. TODO: Add optional install params to setup.py?

"""
import logging
import argparse
import jetstream

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='jetstream draw',
        description=__doc__.replace('``', '"'),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        'path',
        help='Path to a workflow file'
    )

    parser.add_argument(
        'out',
        help='Path to save the image file'
    )

    parser.add_argument(
        '--workflow-format',
        choices=[None, 'yaml', 'json', 'pickle'],
        default=None,
        help='Set the workflow file format instead of using the file '
             'extension.'
    )

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug(args)

    workflow = jetstream.load_workflow(args.path, args.workflow_format)
    workflow.draw(filename=args.out)


if __name__ == '__main__':
    main()
