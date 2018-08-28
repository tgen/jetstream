"""Shared utilities"""
from jetstream.utils.general import *


data_loaders = {
    '.txt': table_to_records,
    '.csv': table_to_records,
    '.mer': table_to_records,
    '.tsv': table_to_records,
    '.json': json_load,
    '.yaml': yaml_load,
    '.yml': yaml_load,
}


def load_data_file(path):
    """Attempts to load a data file from path, raises :ValueError
    if an suitable loader function is not found in data_loaders"""
    for ext, fn in data_loaders.items():
        if path.endswith(ext):
            loader = fn
            break
    else:
        raise ValueError('No loader fn found for {}'.format(path))

    return loader(path)


def loadable_files(directory):
    """Generator yields all files we can load (see data_loaders) """
    for file in os.listdir(directory):
        path = os.path.join(directory, file)
        if os.path.isfile(path) \
                and path.endswith(tuple(data_loaders.keys())):
            yield path


