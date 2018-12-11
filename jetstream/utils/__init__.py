"""Shared utilities"""
from jetstream.utils.general import *


data_loaders = {
    '.txt': load_table,
    '.csv': load_table,
    '.mer': load_table,
    '.tsv': load_table,
    '.json': load_json,
    '.yaml': load_yaml,
    '.yml': load_yaml,
}


def load_data_file(path, filetype=None):
    """Attempts to load a data file from path, raises :ValueError
    if an suitable loader function is not found in data_loaders"""
    if filetype is not None:
        try:
            loader = data_loaders[filetype]
        except KeyError:
            raise ValueError(f'No loader for {filetype}')
    else:
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


