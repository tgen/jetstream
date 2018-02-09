"""Command line utility for generating project reports """
import logging
from jetstream import reports

log = logging.getLogger(__name__)

def arg_parser(subparser):
    parser = subparser.add_parser('report', description=__doc__)
    parser.set_defaults(action=main)

    parser.add_argument('project', nargs='+')

    parser.add_argument('--fast', action='store_true', default=False)

    parser.add_argument('--all', action='store_true', default=False,
                        help='Report on all projects even if they\'re complete.'
                             'Skips complete projects by default.')

    parser.add_argument('--all-jobs', action='store_true', default=False,
                        help='Report on all jobs even if they\'re complete'
                             'Skips complete jobs by default.')

def main(args):
    log.debug(str(args))

    projects = []
    for proj in args.project:
        try:
            p = reports.legacy.Project(proj)
            if p.is_complete and not args.all_jobs:
                log.debug('Skipping complete {}'.format(proj))
                continue
            projects.append(p)
        except FileNotFoundError as err:
            log.critical('Error loading project: {}\n{}'.format(p, err))

    log.debug('Reporting on: {}'.format(str(projects)))
    for p in projects:
        print(reports.legacy.build_plain_text_report(
            p,
            fast=args.fast,
            all_jobs=args.all_jobs
        ))
