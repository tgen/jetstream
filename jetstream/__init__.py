import pkg_resources
from os import environ, listdir, path

_plugin_path = pkg_resources.resource_filename('jetstream', 'plugins/')
environ["JESTREAM_plugin_PATH"] = _plugin_path
plugin_path = _plugin_path
plugins = {f: path.join(_plugin_path, f) for f in listdir(_plugin_path)}

from jetstream.workflow import Workflow
