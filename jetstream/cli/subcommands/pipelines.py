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
from copy import deepcopy
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

    project = jetstream.Project()

    if args.workflow:
        workflow = jetstream.load_workflow(args.workflow)
    else:
        workflow = project.workflow()

    log.info('Project: {} Workflow: {}'.format(project, workflow))

    data = deepcopy(project.config)
    data['project'] = project

    # Load existing project and workflow
    kvarg_data = vars(shared.parse_kvargs(
        args=remaining,
        type_separator=args.kvarg_separator
    ))

    data.update(kvarg_data)

    log.debug('Final template render data:\n{}'.format(data))

    # Render the the new workflow template
    templates = jetstream.render_templates(
        *args.templates,
        data=data,
        search_path=args.search_path
    )

    final_render = '\n'.join(templates)

    log.debug('Final rendered workflow:\n{}'.format(final_render))

    # Combine the old workflow with the new tasks
    if workflow is None:
        # If no workflow is present in the project, just build one
        workflow = jetstream.build_workflow(final_render)
    else:
        # If a workflow is present, we may need to link dependencies in the
        # new workflow to tasks in that workflow, so the build process is
        # different.
        tasks_data = jetstream.utils.yaml_loads(final_render)
        tasks = [jetstream.Task(**data) for data in tasks_data]

        with workflow:
            for t in tasks:
                if not t in workflow:
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

    runner = jetstream.Runner(
        backend=backend,
        max_concurrency=args.max_forks,
        autosave=args.autosave
    )

    runner.start(workflow=workflow, project=project)
    jetstream.save_workflow(workflow, project.workflow_file)

    rc = shared.finalize_run(workflow)
    sys.exit(rc)


if __name__ == '__main__':
    main()
