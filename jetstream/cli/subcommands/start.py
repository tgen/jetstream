"""This module contains the functions that will used to start a new run"""


# TODO add validation to arguments using functions in utils
def validator(value):
    """ Example argument validator function that can be applied with 'type='"""
    if 'condition_is_met':
        return value
    else:
        raise Exception('failed validation')


def arg_parser(subparser):
    parser = subparser.add_parser('start')
    parser.set_defaults(action=main)
    parser.add_argument('config',
                        help='Path to a config file',
                        type=validator)

def main(args):
    print(args)
