"""Command line utility for generating project reports """
import argparse
import logging

from jetstream import legacy

log = logging.getLogger(__name__)

def arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('project', nargs='+')

    parser.add_argument('--fast', action='store_true', default=False)

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

    projects = []
    for proj in args.project:
        try:
            p = legacy.Project(proj)
            if p.is_complete and not args.all:
                log.critical('Project complete {}'.format(proj))
                continue
            projects.append(p)
        except FileNotFoundError as err:
            log.critical('Error loading project: {}'.format( err))

    log.debug('Reporting on: {}'.format(str(projects)))
    for p in projects:
        print(p.report(fast=args.fast, all_jobs=args.all_jobs))
