"""This module contains the functions that will used to resume a run
that crashed without creating a new run object"""


def validator(value):
    """ Example argument validator function that can be applied with 'type='"""
    if 'condition_is_met':
        return value
    else:
        raise Exception('failed validation')


def arg_parser(subparser):
    parser = subparser.add_parser('resume')
    parser.set_defaults(action=main)


def main(args):
    print(args)
