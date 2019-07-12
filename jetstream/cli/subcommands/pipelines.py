"""Run a pipeline.

Pipelines are Jetstream templates that have been documented with version
information and added to a jetstream pipelines directory. This command
allows pipelines to be referenced by name and automatically includes any
pipeline scripts and constants in the run. """
import logging
import os
import jetstream
from jetstream.cli.subcommands import run

log = logging.getLogger('jetstream.cli')
__doc__ = __doc__+ '\n\n' + run.__doc__


def arg_parser(parser):
    # This subcommand is just an extension of the run command, so we add all
    # the arguments from that command to the parser as well.
    run.arg_parser(parser)

    # Hack the existing parser so that it doesn't error if a pipeline
    # name is not given. This allows us to list options instead of getting
    # a parser error immediately.
    for action in parser._actions:
        if 'path' in action.dest:
            action.nargs = '?'

    pipelines = parser.add_argument_group('pipeline options')

    pipelines.add_argument(
        '--searchpath',
        help='override path to the pipelines searchpath',
    )

    pipelines.add_argument(
        '-L', '--list',
        action='store_true',
        help='show a list of all the pipelines installed'
    )


def main(args):
    log.debug(f'{__name__} {args}')

    if args.searchpath:
        jetstream.settings['pipelines']['searchpath'] = args.searchpath

    if args.path:
        # splits the pipeline at @ and takes remaining value as version
        pipeline, *version = args.path.rsplit('@', 1)
        version = next(iter(version), None)

        args.pipeline = pipeline = jetstream.get_pipeline(pipeline, version)
        args.path = os.path.join(pipeline.path, pipeline.info['main'])

        # If bin is declared in the pipeline manifest it should be prepended
        # to the evironment variable PATH
        if pipeline.info.get('bin'):
            bin_path = os.path.join(pipeline.path, pipeline.info['bin'])
            os.environ['PIPELINE_BIN'] = str(bin_path)
            os.environ['PATH'] = f'{bin_path}:{os.environ["PATH"]}'

        run.main(args)
    else:
        if args.list:
            ps = jetstream.list_pipelines()
            ps = '\n'.join(f'{p.name} ({p.version}): {p.path}' for p in ps)
            print(ps)
        else:
            log.error(f'No pipeline name given!')
            raise SystemExit(1)
