"""TODO: update main docstring
"""
import pkg_resources

# I'm adding some package wide variables here in order to cleanup the
# namespace of some subpackages.
plugin_dir = pkg_resources.resource_filename('jetstream', 'plugins/')
plugin_id_pattern = r'(?P<plugin>[^\/]*)\/(?P<path>[^:]*):?(?P<revision>(?<=:)[0-9a-f]{5,40})?$'

from jetstream.workflow import Workflow
from jetstream import formats, plugins
