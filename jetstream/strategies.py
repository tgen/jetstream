import os
import shutil
import tempfile
import subprocess
import hashlib
import logging
import abc

from jetstream import plugins, utils

log = logging.getLogger(__name__)

# TODO probably move the datastores to a different module
def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()


class DataStore(object):
    """ Interface definition for the data store """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def add(self, path):
        raise NotImplementedError

    @abc.abstractmethod
    def cleanup(self):
        raise NotImplementedError


class LocalDataStore(DataStore):
    """ Local data stores rely on a shared filesystem. They initiate a
    temporary directory inside the shared dir. """
    def __init__(self, shared_temp_dir=None):
        self._temp = tempfile.TemporaryDirectory(dir=shared_temp_dir)
        self.path = self._temp.name

    def add(self, source, dest=None):
        if dest is None:
            dest = os.path.join(self._temp.name, source)
        else:
            dest = os.path.join(self._temp.name, dest)

        os.makedirs(os.path.dirname(dest), exist_ok=True)

        m_source = md5sum(source)
        log.debug('Source: {} md5: {}'.format(source, m_source))

        shutil.copy2(source, dest)

        m_dest = md5sum(dest)
        log.debug('Dest: {} md5: {}'.format(dest, m_dest))

        assert m_source == m_dest

    def cleanup(self):
        self._temp.cleanup()


def local(plugin_id):
    """ local strategy executes scripts on the current host. This
     is the default strategy """
    plugin_obj = plugins.get_plugin(plugin, path, revision)

    # These checks should be performed long before we get here
    assert plugin_obj['script'].startswith('#!')
    assert plugin_obj['shell']

    shell = plugin_obj['shell']
    cmd_args = [shell, plugin_obj['_script_path']]

    subprocess.Popen(cmd_args)

