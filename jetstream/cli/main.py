import logging
import argparse
from types import ModuleType

from jetstream.cli import subcommands
from jetstream import __doc__ as readme

log = logging.getLogger(__name__)

def create_parser():
    main_parser = argparse.ArgumentParser(description=readme)

    main_parser.add_argument('--log-level', default='WARNING')


    subparsers = main_parser.add_subparsers(
        title='actions',
        dest='action',
    )

    for _, obj in vars(subcommands).items():
        if isinstance(obj, ModuleType):
            obj.arg_parser(subparsers)

    return main_parser

# TODO This model makes adding subcommands very easy, but it seriously
# slows down the application startup time. Every subcommand module is
# imported for any subcommand used. A solution would be to use parse_known_args
# to gather the args that main knows about, start up logging with those
# settings, then import the subcommand given by args.action and send the
# remaining args to that subcommand.
#
# Current model:
# [rrichholt@it5687:~]$ benchmark.py jetstream
# n Trials: 281
# min: 0.48462734300119337
# max: 0.7841514369938523
# mean: 0.5199922089997532
# median: 0.5213858880015323
# stdev: 0.022948733463308866
#
# With no subcommands imported:
# [rrichholt@it5687:~]$ benchmark.py jetstream
# n Trials: 108
# min: 0.205734956995002
# max: 0.21848879600293003
# mean: 0.2095991600741669
# median: 0.20927944750292227
# stdev: 0.002359137423742718

def main(args=None):
    parser = create_parser()
    args = parser.parse_args(args)

    logging.basicConfig(
        level=getattr(logging, args.log_level)
    )

    if args.action is None:
        parser.print_help()
        print('jetstream -h to see this message')
    else:
        log.debug(str(args))
        args.action(args)
