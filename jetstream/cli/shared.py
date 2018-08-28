import os
import argparse
import jetstream
from jetstream import log


kvarg_types = {
    "default": str,
    "str": str,
    "int": int,
    "float": float,
    "json": jetstream.utils.json_loads,
    "yaml": jetstream.utils.yaml_loads,
    "file": jetstream.load_data_file,
}


def parse_kvargs(args, type_separator=':', types=kvarg_types):
    """Reparses list of arbitrary arguments ``--<type>:<key> <value>``

    This works by building an argument parser configured specifically for the
    arguments present in the list. First any arguments that start with
    ``--`` are selected for creating a new parsing handler. If
    ``type_separator`` is present in the raw argument, it's split to determine
    type and key. The default type for arguments without the separator is
    ``str``. Then an argument is added to the parser with the given type
    and key ().

    After building the parser, it's used to reparse the args and the namespace
    is returned. """
    log.debug('Reparsing kvargs: {}'.format(args))
    parser = argparse.ArgumentParser(add_help=False)

    # TODO another layer of type declaration?
    # --file:json:read_groups ./actually_json.txt
    # or
    # --file-json:read_groups ./actually_json.txt


    for arg in args:
        if arg.startswith('--'):

            if type_separator in arg:
                argtype, _, key = arg.lstrip('-').partition(type_separator)
            else:
                argtype = 'default'
                key = arg.lstrip('-')

            log.debug('Adding parser entry for key: "{}" type: "{}"'.format(
                key, argtype))
            fn = types[argtype]
            
            parser.add_argument(arg, type=fn, dest=key)

    return parser.parse_args(args)
