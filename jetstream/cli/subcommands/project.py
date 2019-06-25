"""View project info or history
"""
import logging
import jetstream

log = logging.getLogger('jetstream.cli')


def arg_parser(parser):
    parser.add_argument(
        '-H', '--history',
        action='store_true'
    )


def project_summary(project):
    wf = project.load_workflow()
    return {
        'index': project.info,
        'paths': vars(project.paths),
        'workflow': {
            'tasks': len(wf),
            'status': wf.summary()
        }
    }


def main(args):
    log.debug(f'{__name__} {args}')

    if args.project:
        if args.history:
            for item in args.project.history_iter():
                print(jetstream.utils.yaml_dumps(item))
        else:
            summary = project_summary(args.project)
            print(jetstream.utils.yaml_dumps(summary))

    else:
        err = 'No project given! Must be run inside a project, or use ' \
              '-p/--project argument'
        raise ValueError(err)
