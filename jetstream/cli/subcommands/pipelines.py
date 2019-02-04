"""Run a pipeline

Pipelines are Jetstream templates that have been documented with version
information and added to the jetstream pipelines directory. This command
allows pipelines to be referenced by name and automatically includes the
pipeline scripts and constants in the run."""
import os
import logging
import jetstream
from jetstream.cli.subcommands.run import main as run_main
from jetstream.cli.subcommands.run import arg_parser as run_arg_parser

log = logging.getLogger(__name__)


def arg_parser(parser):
    # This subcommand is just an extension of the run command, so we add all
    # the arguments from that command to the parser as well.
    run_arg_parser(parser)

    parser.add_argument(
        '--pipelines-home',
        help='Override path to the pipelines.',
    )

    parser.add_argument(
        '--ls',
        action='store_true',
        help='Show details about installed pipelines'
    )


def main(args=None):
    log.debug(f'{__name__} {args}')

    # Make sure pipelines are configured
    if args.pipelines_home:
        jetstream.settings['pipelines']['home'] = args.pipelines_home

    pipelines_home = jetstream.settings['pipelines']['home'].get()
    if not pipelines_home or not os.path.isdir(pipelines_home):
        err = 'Pipelines are not configured. Check the application settings.'
        raise ValueError(err)

    # Ls just prints info about current pipelines
    if args.ls:
        installed = jetstream.pipelines.ls(pipelines_home)
        return print(jetstream.utils.yaml_dumps(installed))

    # If the pipeline has extra constants, add them to the settings constants
    pipeline = jetstream.pipelines.lookup(args.path)
    args.path = os.path.join(pipeline.path, pipeline.manifest['main'])
    jetstream.settings['constants'].set(pipeline.manifest.get('constants'))

    # If the pipeline has a bin directory, prepend the env PATH
    if 'bin' in pipeline.manifest:
        bin_path = os.path.join(pipeline.path, pipeline.manifest['bin'])
        os.environ['PIPELINE_BIN'] = str(bin_path)
        os.environ['PATH'] = f'{bin_path}:{os.environ["PATH"]}'

    run_main(args)


if __name__ == '__main__':
    main()
