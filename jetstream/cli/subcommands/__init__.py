# This package contains all  of the argument parsers and main functions
# for Jetstream commands. When adding to this package, please follow the
# template below to give the subcommands some consistent behaviors. They
# are dynamically imported in the jetstream.cli.__init__ module and will need
# to be added to the _subcommands dict there to be found.
#
#
# """Short help text for subcommand
# Any additional lines are only shown for the command help
# """
# import logging
#
# log = logging.getLogger(__name__)
#
# def add_arguments(parser):
#     # Add any command arguments to the parser here
#
# def main(args):
#     log.debug(f'{__name__} {args}')
#     # Here is the main function for this command, it receives the argparse
#     # Namespace object as the only argument. Args will always have an
#     # additional attribute "kvargs" which contains any extra workflow data
#     # passed in as command-line args.

# Notes:
# run, render, build, and pipelines all share some common options which are
# included as a separate module. The