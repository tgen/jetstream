import os
import argparse
import jetstream
from jetstream import log


kvarg_types = {
    "default": str,
    "str": str,
    "int": int,
    "float": float,
    "bool": jetstream.utils.to_bool,
    "json": jetstream.utils.json_loads,
    "yaml": jetstream.utils.yaml_loads,
    "file": jetstream.load_data_file,
}


def parse_kvargs(args, separator=None):
    """Reparses list of arbitrary arguments ``--<type>:<key> <value>``

    This works by building an argument parser configured specifically for the
    arguments present in the list. First any arguments that start with
    ``--`` are selected for creating a new parsing handler. If
    ``type_separator`` is present in the raw argument, it's split to determine
    type and key. Then an argument is added to the parser with the given type
    and key.

    After building the parser, it's used to reparse the args list and the
    namespace is returned as a dictionary. """
    log.debug('Reparsing kvargs: {}'.format(args))
    parser = argparse.ArgumentParser(add_help=False)

    # TODO another layer of type declaration?
    # --file:json:read_groups ./actually_json.txt
    # or
    # --file-json:read_groups ./actually_json.txt

    if separator is None:
        separator = ':'

    for arg in args:
        if arg.startswith('--'):
            argtype, _, key = arg.lstrip('-').partition(separator)

            if argtype and key:
                if argtype not in kvarg_types:
                    raise ValueError(
                        f'Variable: "{arg}" unknown type: "{argtype}"'
                    )

                log.debug(f'Adding parser for key: "{key}" type: "{argtype}"')
                fn = kvarg_types[argtype]
                parser.add_argument(arg, type=fn, dest=key)

    return vars(parser.parse_args(args))


def set_project(path=None):
    """If path is set, chdir to the given project and load the object.
    Otherwise check if the cwd is a project and return or return None."""
    if path:
        os.chdir(path)
        project = jetstream.Project()
    else:
        try:
            project = jetstream.Project()
        except jetstream.NotAProject:
            project = None

    return project


def load_variables(path):
    """Most variables files should load with the yaml parser, but tabs in
    a json file might raise an error, so both are tried."""
    try:
        return jetstream.utils.yaml_load(path)
    except jetstream.utils.yaml.YAMLError:
        return jetstream.utils.json_load(path)
