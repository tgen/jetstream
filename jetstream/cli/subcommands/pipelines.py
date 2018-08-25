"""Run a Jetstream workflow inside a project

Template variable data is usually saved as files in ``<project>/config``, but 
command arguments can also be used to pass variable data to templates. All 
arguments remaining after parsing the command line (arguments that are not 
listed in the help section) will be treated as template variable data (kvargs):

Template variable arguments should follow the syntax: ``--<key> <value>``. 
The key must start with two hyphens and the value is the following argument. The 
variable type can be explicitly set with the syntax ``--<type>:<key> <value>``.
Variables with no type declared will be loaded as strings.

If the variable type is "file" the value will be ``jetstream.data_loaders``,
which handles files according to their extension. All other types will
evaluated by the appropriate type function.

"""
import sys
import jetstream
from jetstream import log
from jetstream.cli import shared
from jetstream.cli.subcommands.run import arg_parser


def generate_workflow(templates, data, search_path, render_only, build_only):
    log.debug('Template render data: {}'.format(data))

    templates = jetstream.render_templates(
        *templates,
        data=data,
        search_path=search_path
    )

    if render_only:
        for t in templates:
            print(t)
        sys.exit(0)

    workflow = jetstream.build_workflow('\n'.join(templates))

    if build_only:
        print(workflow.to_yaml())
        sys.exit(0)

    return workflow


def main(args=None):
    parser = arg_parser()
    parser.prog = 'jetstream pipelines'
    args, remaining = parser.parse_known_args(args)
    log.debug(args)

    project = jetstream.Project()
    workflow = project.workflow()

    log.info('Project: {} Workflow: {}'.format(project, workflow))

    data = vars(shared.parse_kvargs(
        args=remaining,
        type_separator=args.kvarg_separator
    ))

    if 'project' not in data:
        data['project'] = project


    new_workflow = generate_workflow(
        templates=args.templates,
        data=data,
        search_path=args.search_path,
        render_only=args.render_only,
        build_only=args.build_only
    )

    if workflow is None:
        workflow = new_workflow
    else:
        workflow.compose(new_workflow)
    
    if args.method == 'retry':
        workflow.retry()
    else:
        workflow.resume()

    log.info('Workflow data after composition: {}'.format(workflow))

    if args.backend == 'slurm':
        backend = jetstream.SlurmBackend()
    else:
        backend = jetstream.LocalBackend()

    runner = jetstream.AsyncRunner(
        backend=backend,
        max_forks=args.max_forks,
        autosave=args.autosave
    )

    rc = runner.start(workflow=workflow, project=project)
    runner.close()

    jetstream.save_workflow(workflow, project.workflow_file)
    sys.exit(rc)


if __name__ == '__main__':
    main()
