#!/usr/bin/env python3
import logging
import sys

from pydot import graph_from_dot_file

import jetstream.workflow.async_runner

# TODO These should be written with a unittest framework

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    log.critical('Logging started')

    graph = graph_from_dot_file(sys.argv[1])[0]
    wf = jetstream.Workflow.from_pydot(graph)
    jetstream.workflow.async_runner.run(wf, debug=False)
