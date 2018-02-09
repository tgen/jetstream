"""Command line utility for managing plugins """
import logging
from os import getcwd, listdir
from jetstream import reports

log = logging.getLogger(__name__)

def arg_parser(subparser):
    parser = subparser.add_parser('report', description=__doc__)
    parser.set_defaults(action=main)

    parser.add_argument('project', nargs='*', default=list(find_projects()))

    parser.add_argument('--fast', action='store_true', default=False)

    parser.add_argument('--all-jobs', action='store_true', default=False)


def find_projects():
    """ Yields incomplete projects in the cwd """
    for d in listdir(getcwd()):
        try:
            p = reports.legacy.Project(d)
            if not p.is_complete:
                yield d
        except FileNotFoundError:
            pass
    raise StopIteration


def main(args):
    log.debug(str(args))

    projects = []
    for p in args.project:
        try:
            projects.append(reports.legacy.Project(p))
        except FileNotFoundError as err:
            log.critical('Error loading project: {}\n{}'.format(p, err))

    log.debug('Reporting on: {}'.format(str(projects)))
    for p in projects:
        print(reports.legacy.build_plain_text_report(
            p,
            fast=args.fast,
            all_jobs=args.all_jobs
        ))
