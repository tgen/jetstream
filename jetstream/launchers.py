import os
import shutil
import tempfile
import subprocess
import hashlib
import logging
import random
from datetime import datetime
from jetstream.workflow import Result
from jetstream.utils import fingerprint
import abc

log = logging.getLogger(__name__)


def dry(plugin):
    """ This dry-run strategy can serve for testing and also as a model
    for other, more complex, strategies. I'm going to add some notes here
    for development of other strategies.

    # Plugins

    Plugins exist in two states:

    1) the external representation, yaml serialization that contains all of the
    data that will be loaded and parsed by the plugins module.

    2) the internal representation, a dictionary as parsed by the yaml reader

    These functions (strategies) should handle the latter. And they should
    return a record of the plugin execution (TBD).

    # Scripts

    The core of a plugin is the script. At the present time, a script is the
    only requirement for a plugin. But, plugins will eventually exist which
    require diverse execution paradigms. It is very important that the external
    representation (1) of a plugin can fully accommodate these paradigms. If we
    don't get this right from the start, we'll constantly need to modify the
    plugin format, and they become entangled with the application development.

    So.. what should the constraints be for a script? They need to be platform-
    neutral instruction sets. Here are the current options I see:

    Bash works well for this because it is available on virtually any system.
    But, we don't want to limit it to just bash, because there are a lot of
    modules in this codebase that enable python scripts to work with jetstream
    projects easier.

    Another attractive option is allowing shebangs to determine the program used
    to execute the script. But, consider this:

    The shebang line has never been specified as part of POSIX, SUS, LSB or any
    other specification. AFAIK, it hasn't even been properly documented.

    From Jörg W Mittag on Stack Overflow
    (https://stackoverflow.com/a/4304187/3924113):

    "There is a rough consensus about what it does: take everything between the
    ! and the \n and exec it. The assumption is that everything between the !
    and the \n is a full absolute path to the interpreter. There is no consensus
     about what happens if it contains whitespace.

    1. Some operating systems simply treat the entire thing as the path. After
    all, in most operating systems, whitespace or dashes are legal in a path.

    2. Some operating systems split at whitespace and treat the first part as
    the path to the interpreter and the rest as individual arguments.

    3. Some operating systems split at the first whitespace and treat the front
    part as the path to the interpeter and the rest as a single argument (which
    is what you are seeing [if you try this #!/usr/bin/gawk --re-interval -f]).

    4. Some even don't support shebang lines at all.

    Thankfully, 1. and 4. seem to have died out, but 3. is pretty widespread, so
    you simply cannot rely on being able to pass more than one argument.

    And since the location of commands is also not specified in POSIX or SUS,
    you generally use up that single argument by passing the executable's name
    to env so that it can determine the executable's location; e.g.:

    #!/usr/bin/env gawk
    [Obviously, this still assumes a particular path for env, but there are only
    very few systems where it lives in /bin, so this is generally safe. The
    location of env is a lot more standardized than the location of gawk or even
    worse something like python or ruby or spidermonkey.]

    Which means that you cannot actually use any arguments at all."

    This means that there is no guarantee that a script with a shebang behaves
    the same from system to system. Although, I am not convinced that is a
    defeating problem (the code inside the scripts themselves is more prone to
    system variability than the shebang), it should be recognized.

    It is most likely that some attribute will be added to plugins to control
    the interpreter to be used. In this event, plugin repos will probably exist
    under their own systems of control, where users submit plugins intended
    to be run under a certain environment. The administrator of the environment
    would then be responsible for validating plugins and make a requirement
    for the interpreters that can be used.

    """
    log.critical('Starting dry run for {}'.format(plugin['id']))

    result = {
        'id': plugin['id'],
        'rc': 0,
        'started': str(datetime.now()),
        'completed': str(datetime.now()),
        'fingerprint': fingerprint(),
        'log': 'Dry run'
    }

    return result


def default(plugin):
    """Launches plugin by guessing the interpreter to use from shebang. This
    currently supports Python2/3 and Bash. """
    log.critical('Starting plugin {}'.format(plugin['id']))

    shebang = plugin['script'].splitlines()[0]

    if not shebang.startswith('#!'):
        log.warning('Unable to parse shebang in {}'.format(plugin['id']))
        shell_cmd = ['bash']
    elif 'python3' in shebang:
        shell_cmd = ['python3']
    elif 'python' in shebang:
        shell_cmd = ['python']
    else:
        shell_cmd = ['bash']

    p = subprocess.Popen(
        shell_cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    stdout, _ = p.communicate(input=plugin['script'].encode())

    result = Result(
        plugin_id=plugin['id'],
        return_code=p.returncode,
        log=stdout.decode()
    )

    return result


def beta_stageinout(plugin):
    """ This is an early example of how we might stage project data in temporary
    data stores for remote execution. """
    log.critical('Starting run for {}'.format(plugin['id']))

    # Setup phase
    work_dir = LocalFsDataStore()
    log.critical('Setup workspace in {}'.format(work_dir.path))

    stagein = plugin.get('stage_in')
    if stagein:
        subprocess.check_call(['rsync', '--progress', '--exclude',
                               '.jetstream', '-rtu', '.', work_dir.path])

    try:
        work_dir.add('project.json')
    except FileExistsError:
        pass

    # Execution phase
    log.debug('Launching {}'.format(plugin['id']))
    rnd = random.randrange(30, 45)
    cmd = 'echo $(date) $(pwd) "{}" | tee logs{}.bam && sleep {}'.format(
        plugin['id'], rnd, rnd)
    dry_script = ['bash', '-vx', '-c', cmd]
    p = subprocess.Popen(
        dry_script,
        cwd=work_dir.path,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE
    )
    p.wait()

    # Teardown phase
    stageout = plugin.get('stage_out')
    log.critical('Copying {} back from {}'.format(stageout, work_dir.path))

    if stageout:
        subprocess.check_call(['rsync', '--progress', '--exclude',
                               'project.json', '-rtu', work_dir.path + '/', '.'])

    # Generate results
    # TODO The structure of this object is determined by Workflow.__send__()
    result = {
        'plugin_id': plugin['id'],
        'stdout': p.stdout.read().decode(),
        'stderr': p.stderr.read().decode(),
        'rc': p.returncode
    }

    return result


def md5sum(filename, blocksize=65536):
    """Returns a md5 checksum for a given file. This processes the data
    in chunks set by 'blocksize'. """
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()


def safe_copy(source, dest):
    """Copies a source to dest and verifies md5 checksum integrity"""
    os.makedirs(os.path.dirname(dest), exist_ok=True)

    log.debug('Copy: {} -> {}'.format(source, dest))

    md5_source = md5sum(source)
    log.debug('Source: {} md5: {}'.format(source, md5_source))

    shutil.copy2(source, dest)

    md5_dest = md5sum(dest)
    log.debug('Dest: {} md5: {}'.format(dest, md5_dest))

    assert md5_source == md5_dest


class DataStore(object):
    """ Class definition for the data store """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def add(self, source, dest=None):
        raise NotImplementedError

    @abc.abstractmethod
    def cleanup(self):
        raise NotImplementedError


class LocalFsDataStore(DataStore):
    def __init__(self, *args, **kwargs):
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
    def __init__(self, shared_temp_dir=None):
        self._temp = tempfile.TemporaryDirectory(dir=shared_temp_dir)
        self.path = self._temp.name

    def add(self, source, dest=None):
        if dest is None:
            dest = os.path.join(self._temp.name, source)
        else:
            dest = os.path.join(self._temp.name, dest)

        safe_copy(source, dest)

    def cleanup(self):
        self._temp.cleanup()
