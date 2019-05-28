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
        'config': args.project.get_config(),
        'info': args.project.get_info(),
        'workflow': str(args.project.get_workflow())
    }
    print(jetstream.utils.yaml_dumps(p))
