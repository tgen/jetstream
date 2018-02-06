"""TODO: update main docstring


"""
import pkg_resources
import logging

log = logging.getLogger(__name__)

from . import utils

# I'm adding some package wide variables here in order to cleanup the
# namespace of some subpackages.
PLUGIN_DIR = pkg_resources.resource_filename('jetstream', 'plugins/')
PLUGIN_ID_PATTERN = r'(?P<plugin>[^\/]*)\/(?P<path>[^:]*):?(?P<revision>(?<=:)[0-9a-f]{5,40})?$'


# TODO: should project functions walk up the directory tree like git?
# see this https://gist.github.com/zdavkeos/1098474


from . import formats, plugins, config, project
from .workflow import Workflow
