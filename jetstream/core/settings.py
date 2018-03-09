"""Contains global application settings"""
from ruamel import yaml
from os import environ

DEFAULT = {
    'RUN_DATA_DIR': '.jetstream',
    'RUN_DATA_PREFIX': 'js'
}

profile_path = environ.get('JETSTREAM_PROFILE')

if profile_path is not None:
    with open(profile_path, 'r') as fp:
        profile = yaml.load(stream=fp)
else:
    profile = DEFAULT
