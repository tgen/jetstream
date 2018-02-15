import json
import logging

from jetstream.utils import yaml
from .validator import check

log = logging.getLogger(__name__)


def deserialize(data, format='yaml'):
    if format == 'yaml':
        config = yaml.load(data, Loader=yaml.Loader)
    elif format == 'json':
        config = json.loads(data)
    else:
        config = format(data)

    return config


def read(path, *args, **kwargs):
    """ Read a config file from a path. Path can be "yaml", "json", or a
     function that takes data as its first argument and returns a config. """
    with open(path, 'r') as fp:
        data = fp.read()

    return deserialize(data, *args, **kwargs)


def serialize(config, format='yaml'):
    if format == 'yaml':
        data = yaml.dump(config, default_flow_style=False)
    elif format == 'json':
        data = json.dumps(config, indent=4)
    else:
         data = format(config)
    return data


def write(config, path, *args, **kwargs):
    with open(path, 'w') as fp:
        fp.write(serialize(config, *args, **kwargs))


def validate(config):
    """ Validate the key-value pair against the schema"""
    return check(config['meta'])