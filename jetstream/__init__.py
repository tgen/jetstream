"""Jetstream is a collection of tools for automating workflows at TGen."""
import pkg_resources

# I'm adding some package wide variables here in order to cleanup the
# namespace of some subpackages.
PLUGIN_DIR = pkg_resources.resource_filename('jetstream', 'plugins/')
PLUGIN_ID_PATTERN = r'(?P<plugin>[^\/]*)\/(?P<path>[^:]*):?(?P<revision>(?<=:)[0-9a-f]{5,40})?$'


from . import formats, plugins, config, project, launch, utils
from .workflow import Workflow

# TODO: Move reports.legacy.Project into jetstream.project

# TODO: should project functions walk up the directory tree like git?
# see this https://gist.github.com/zdavkeos/1098474

