"""Note: This command is primarily intended for debugging template issues. To
render, build, and run a workflow in one step use "jetstream_pipelines". This
command inherits its arguments from jetstream_pipelines. """
import logging
import jetstream
import jetstream_pipelines
from jetstream_pipelines.main import arg_parser

log = logging.getLogger(__name__)

# Argument parser should be nearly identical to jetstream_pipelines.main

def main(args=None):
    parser = arg_parser()
    parser.epilog = __doc__
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    p = jetstream.Project()
    env = jetstream_pipelines.env()
