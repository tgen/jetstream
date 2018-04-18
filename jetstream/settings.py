"""Contains global application settings"""
import yaml
import logging
from os import environ, path

log = logging.getLogger(__name__)

home = path.expanduser('~')

profile_default_path = path.join(home, '.jetstream.yaml')
profile_path = environ.get('JETSTREAM_PROFILE', profile_default_path)

defaults = {
    'JETSTREAM_RELEASES_DIR': path.join(home, '.jetstream_pipelines')
}

profile = defaults.copy()

try:
    if path.exists(profile_path):
        with open(profile_path, 'r') as fp:
            profile.update(yaml.load(fp.read()))
    environ.update(profile)
except Exception as e:
    log.critical('Error loading profile {}: {}'.format(profile_path, e))
    profile = defaults
    environ.update(profile)
