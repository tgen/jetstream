"""Interact with jetstream projects
View tasks, run history, or project data.
"""
import logging
import jetstream

log = logging.getLogger(__name__)


def arg_parser(parser):
    pass


def main(args):
    log.debug(f'{__name__} {args}')
    p = {
        'config': args.project.config,
        'info': args.project.info,
        'workflow': str(args.project.workflow)
    }
    print(jetstream.utils.yaml_dumps(p))
