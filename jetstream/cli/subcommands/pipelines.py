"""Run a pipeline.

Pipelines are Jetstream templates that have been documented with version
information and added to a Jetstream pipelines directory. This command
allows pipelines to be referenced by name and automatically includes any
pipeline scripts and constants in the run.

Pipeline names may also include a version number <name>@<version>.

For complete option listing see "jetstream run" """
import logging
import sys
import errno
import jetstream
from jetstream.cli.subcommands import run_common_options, run

log = logging.getLogger('jetstream.cli')


def add_arguments(parser):
    # # This subcommand is just an extension of the run command, so we add all
    # # the arguments from that command to the parser as well.
    pipelines = parser.add_argument_group('pipelines command options')

    pipelines.add_argument(
        'name',
        nargs='?',  # If not given, all pipelines will be listed
        help='pipeline name (optionally with a version number)'
    )

    pipelines.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='show detailed description of the pipeline'
    )

    run_common_options.add_arguments(parser)


def main(args):
    log.debug(f'{__name__} {args}')
    if args.search_path:
        searchpath = args.search_path
    else:
        searchpath = jetstream.settings['pipelines']['searchpath'].get()
        searchpath = searchpath.split(':')

    # load and describe pipeline
    if args.verbose and args.pipeline:
        # Direct pipeline path given, just load and show details
        ctx = args.pipeline.get_context()
        try:
            jetstream.utils.dump_yaml(ctx, sys.stdout)
        except IOError as e:
            if e.errno == errno.EPIPE:
                pass
    elif args.verbose and args.name:
        # Find a pipeline and show details
        pipeline, *version = args.name.rsplit('@', 1)
        version = next(iter(version), None)
        args.pipeline = jetstream.get_pipeline(pipeline, version, searchpath=searchpath)
        ctx = args.pipeline.get_context()
        try:
            jetstream.utils.dump_yaml(ctx, sys.stdout)
        except IOError as e:
            if e.errno == errno.EPIPE:
                pass
    elif args.pipeline:
        # Direct pipeline path give, just load and run
        run.main(args)
    elif args.name:
        # Find and run a pipeline
        pipeline, *version = args.name.rsplit('@', 1)
        version = next(iter(version), None)
        args.pipeline = jetstream.get_pipeline(pipeline, version, searchpath=searchpath)
        run.main(args)
    else:
        # List all pipelines
        all_pipelines = jetstream.list_pipelines(*searchpath)
        output = '\n'.join(f'{p.name} ({p.version}): {p.path}' for p in all_pipelines)
        try:
            print(output)
        except IOError as e:
            if e.errno == errno.EPIPE:
                pass
