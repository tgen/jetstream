#!/usr/bin/env python3

""" Example running a jetstream workflow that has been saved in dot format

This script loads a workflow from the path given as the first argument, and
runs it with the synchronous workflow runner.

"""
import sys
import logging
from jetstream import workflow

# TODO These should be written with a unittest framework

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    log.critical('Logging started')

    wf = workflow.load_pydot(sys.argv[1])
    workflow.runner.run(wf, debug=False)

    log.critical('All done!')
