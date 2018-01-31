import os
import shutil
import tempfile
import subprocess
import hashlib
import logging
import time
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
    """ Class definition for the data store """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def setup(self):
        raise NotImplementedError

    @abc.abstractmethod
    def add(self, source, dest=None):
        raise NotImplementedError

    @abc.abstractmethod
    def cleanup(self):
        raise NotImplementedError


class LocalFsDataStore(DataStore):
    def setup(self, *args, **kwargs):
        self._temp = tempfile.TemporaryDirectory()
        self.path = self._temp.name

    def add(self, source, dest=None):
        if dest is None:
            dest = os.path.join(self.path, source)
        else:
            dest = os.path.join(self.path, dest)

        os.makedirs(os.path.dirname(dest), exist_ok=True)

        log.debug('Symlink: {} -> {}'.format(source, dest))
        os.symlink(source, dest)
        # TODO This is unix only, can we make this platform independent?

    def cleanup(self):
        self._temp.cleanup()


class SharedFsDataStore(DataStore):
    """Shared filesystem data stores rely on a filesystem being shared between
    the runner host and the strategy hosts. They initiate a temporary directory
    on the shared filesystem. Files added to the store will by copied to the
    temp directory and verified. """
    def setup(self, shared_temp_dir=None):
        self._temp = tempfile.TemporaryDirectory(dir=shared_temp_dir)
        self.path = self._temp.name

    def add(self, source, dest=None):
        if dest is None:
            dest = os.path.join(self._temp.name, source)
        else:
            dest = os.path.join(self._temp.name, dest)

        os.makedirs(os.path.dirname(dest), exist_ok=True)

        log.debug('Copy: {} -> {}'.format(source, dest))

        md5_source = md5sum(source)
        log.debug('Source: {} md5: {}'.format(source, md5_source))

        shutil.copy2(source, dest)

        md5_dest = md5sum(dest)
        log.debug('Dest: {} md5: {}'.format(dest, md5_dest))

        assert md5_source == md5_dest
        # TODO how do we want to handle this error? I have no idea
        # how often it might happen, if at all. If it's extremely
        # rare, then exception handling may not be necessary

    def cleanup(self):
        self._temp.cleanup()


class Local(object):
    def __init__(self, project):
        self.project = project

    def setup(self):
        self.store = LocalFsDataStore()


def local(plugin_id):
    """ local strategy executes scripts on the current host. This
     is the default strategy """
    plugin_obj = plugins.get_plugin(plugin_id)

    # These checks should be performed long before we get here
    assert plugin_obj['script'].startswith('#!')
    assert plugin_obj['shell']

    shell = plugin_obj['shell']
    cmd_args = [shell, plugin_obj['_script_path']]

    subprocess.Popen(cmd_args)


# TODO allow different execution strategies here
# they will essentially be a separate thread, but they can return results

# TODO allow asynchronous strategies
# These will be harder to code, but performance-wise they're better if that
# becomes an issue


def dummy_strategy(plugin):
    print('strategy running with: {}'.format(plugin))


    # Plugin component setup
    # Setup data store


    # Plugin component teardown
    # Cleanup data store

    time.sleep(30)

    record = {
        "pid": 'Threaded strategy',
        "returncode": 0,
        "stdout": 'strategy done {}'.format(plugin),
        "stderr": None,
    }

    return record
