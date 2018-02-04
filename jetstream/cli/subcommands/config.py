""" This module contains the cli interface code for the config file utility.
"""

# TODO:
# Allow cli access to config file validator functions
# Allow cli access to config file converter functions



def validator(value):
    """ Example argument validator function that can be applied with 'type='"""
    if 'condition_is_met':
        return value
    else:
        raise Exception('failed validation')


def arg_parser(subparser):
    parser = subparser.add_parser('config')
    parser.set_defaults(action=main)
    parser.add_argument('example', help='Change me!', type=validator)


def main(args):
    print(args)
