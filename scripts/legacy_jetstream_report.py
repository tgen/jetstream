#!/usr/bin/env python
"""Generate project reports.

This tool only works with Pegasus/Medusa projects."""
import argparse
import logging
from jetstream.legacy.project import Project

log = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('project', nargs='+',
                        help='Path to a Pegasus/Medusa project. Multiple '
                             'projects can be given.')

    parser.add_argument('--fast', action='store_true', default=False,
                        help='Only reports whether or not the project is'
                             'complete. This skips the job status reporting.')

    parser.add_argument('--all', action='store_true', default=False,
                        help='Report on all projects even if they\'re complete.'
                             ' Skips complete projects by default.')

    parser.add_argument('--all-jobs', action='store_true', default=False,
                        help='Report on all jobs even if they\'re complete. '
                             'Skips complete jobs by default.')
    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    log.debug('{}: {}'.format(__name__, args))

    logging.basicConfig(
        level=logging.INFO
    )

    projects = []
    for proj in args.project:
        try:
            p = Project(proj)
            if p.is_complete and not args.all:
                log.info('Project complete {}'.format(proj))
                continue
            projects.append(p)
        except FileNotFoundError as err:
            log.info('Error loading project: {}'.format(err))

    log.debug('Reporting on: {}'.format(str(projects)))
    for p in projects:
        print(p.report(fast=args.fast, all_jobs=args.all_jobs))

