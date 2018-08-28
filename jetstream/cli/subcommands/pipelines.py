"""Run a Jetstream workflow inside a project.

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
from jetstream.backends import LocalBackend, SlurmBackend
from jetstream.cli.subcommands.run import arg_parser as run_arg_parser


def arg_parser():
    parser = run_arg_parser()
    parser.prog = 'jetstream pipelines'
    parser.description = __doc__
    return parser


def main(args=None):
    parser = arg_parser()
    args, remaining = parser.parse_known_args(args)
    log.debug(args)

    # Load existing project and workflow
    variable_data = vars(shared.parse_kvargs(
        args=remaining,
        type_separator=args.kvarg_separator
    ))

    project = jetstream.Project()
    workflow = project.workflow()

    if 'project' not in variable_data:
        variable_data['project'] = project

    log.info('Project: {} Workflow: {}'.format(project, workflow))

    # Render the the new workflow template
    templates = jetstream.render_templates(
        *args.templates,
        data=variable_data,
        search_path=args.search_path
    )

    # Combine the old workflow with the new tasks
    if workflow is None:
        # If no workflow is present in the project, just build one
        workflow = jetstream.build_workflow('\n'.join(templates))
    else:
        # If a workflow is present, we may need to link dependencies in the
        # new workflow to tasks in that workflow, so the build process is
        # different.
        rendered_templates = jetstream.render_template('\n'.join(templates))
        tasks_data = jetstream.utils.yaml_loads(rendered_templates)
        tasks = [jetstream.Task(**data) for data in tasks_data]

        with workflow:
            for t in tasks:
                workflow.add_task(t)

    # Reset the workflow tasks according to method
    if args.method == 'retry':
        workflow.retry()
    else:
        workflow.resume()

    log.info('Workflow data after composition: {}'.format(workflow))

    if args.backend == 'slurm':
        backend = SlurmBackend()
    else:
        backend = LocalBackend()

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
