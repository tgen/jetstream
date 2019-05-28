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

    pipelines = parser.add_argument_group('pipeline options')

    pipelines.add_argument(
        '--pipelines-home',
        help='Override path to the pipelines home',
    )

    pipelines.add_argument(
        '--ls',
        action='store_true',
        help='Show details about installed pipelines'
    )

# TODO Arg parser upgrades?
# All these share an arg parser:
# jetstream pipelines name
# jetstream build name / jetstream run name / jetstream render name
# the value "name" should be optional positional, then:
# jetstream build <args> == jetstream run --build-only <args>
# jetstream render <args> == jetstream run --render-only <args>
# jetstream pipelines name == jetstream run -t pipeline/main.jst <args>
# jetstream pipelines

def main(args=None):
    log.debug(f'{__name__} {args}')

    # Make sure pipelines are configured
    if args.pipelines_home:
        jetstream.settings['pipelines']['home'] = args.pipelines_home

    pipelines_home = jetstream.settings['pipelines']['home'].get(str)

    # ls just prints info about current pipelines
    if args.ls:
        installed = jetstream.pipelines.ls(pipelines_home)
        print(jetstream.utils.yaml_dumps(installed))
    else:
        args.pipeline = pipeline = jetstream.pipelines.lookup(args.path)
        args.path = os.path.join(pipeline.path, pipeline.manifest['main'])

        # If the pipeline has a bin directory, prepend the env PATH
        if 'bin' in pipeline.manifest:
            bin_path = os.path.join(pipeline.path, pipeline.manifest['bin'])
            os.environ['PIPELINE_BIN'] = str(bin_path)
            os.environ['PATH'] = f'{bin_path}:{os.environ["PATH"]}'

        run_main(args)


if __name__ == '__main__':
    main()
