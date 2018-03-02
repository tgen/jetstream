"""Jetstream is a collection of tools for automating workflows at TGen."""
import pkg_resources

__version__ = '0.1.0a1'
__author__ = 'Ryan Richholt'
#TODO finish these attributes

PLUGIN_DIR = pkg_resources.resource_filename('jetstream', 'plugins/')
PLUGIN_ID_PATTERN = r'(?P<plugin>[^\/]*)\/(?P<path>[^:]*):?(?P<revision>(?<=:)[0-9a-f]{5,40})?$'

from .settings import profile, profile_path
