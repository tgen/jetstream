"""TODO: update main docstring
"""
import pkg_resources

# I'm adding some package wide variables here in order to cleanup the
# namespace of some subpackages.
PLUGIN_DIR = pkg_resources.resource_filename('jetstream', 'plugins/')
PLUGIN_ID_PATTERN = r'(?P<plugin>[^\/]*)\/(?P<path>[^:]*):?(?P<revision>(?<=:)[0-9a-f]{5,40})?$'

from . import formats, plugins, config, utils
from .workflow import Workflow
