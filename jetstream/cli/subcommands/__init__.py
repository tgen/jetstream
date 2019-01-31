""""""

# Subcommands module contains all the argument parsers for Jetstream
# commands. When adding to this package, please follow the template
# below to give the subcommands some consistent behaviors:
#
# .. code-block:: python
#
# """Short help text for subcommand
# Any additional lines are only shown for the command help
# """
# import logging
#
# log = logging.getLogger(__name__)
#
# def arg_parser(parser):
#     # Add any command arguments to the parser here
#
#
# def main(args):
#     log.debug(f'{__name__} {args}')
#     # Here is the main function for this command, it receives the argparse
#     # Namespace object as the only argument. Args will always have an
#     # additional attribute "kvargs" which contains any extra workflow data
#     # passed in as command-line args.
#
