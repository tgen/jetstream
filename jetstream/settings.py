"""Contains global application settings"""
import logging
from os import environ, path
from ruamel import yaml

log = logging.getLogger(__name__)

home = path.expanduser('~')

profile_default_path = path.join(home, '.jetstream.yaml')
profile_path = environ.get('JETSTREAM_PROFILE', profile_default_path)

defaults = {}

profile = defaults.copy()

try:
    if path.exists(profile_path):
        with open(profile_path, 'r') as fp:
            profile.update(yaml.load(fp.read()))
        log.debug('Loaded profile: {}'.format(profile_path))
    else:
        log.debug('Loaded default profile')
except Exception as e:
    log.critical('Error loading profile {}: {}'.format(profile_path, e))


environ.update(profile)
