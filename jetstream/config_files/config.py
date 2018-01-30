import json
from ruamel import yaml
import logging

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


def write(path, config, *args, **kwargs):
    with open(path, 'w') as fp:
        fp.write(serialize(config, *args, **kwargs))


def validate(config):
    """ Validate the key-value pair against the schema"""
    fails = []
    for k, v in config['meta'].items():
        if check(k, v):
            pass
        else:
            fails.append((k, v))

    if fails:
        for k, v in fails:
            log.warning('Validator failed for "{}: {}"!'.format(k, v))
        raise ValueError(str(fails))