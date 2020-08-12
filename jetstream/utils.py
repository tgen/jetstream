"""Shared utilities"""
import builtins
import confuse
import csv
import fnmatch
import gzip
import importlib
import io
import json
import logging
import os
import subprocess
import sys
import yaml
from collections.abc import Sequence, Mapping
from datetime import datetime
from getpass import getuser
from multiprocessing import cpu_count
from socket import gethostname
from uuid import getnode, uuid4

import jetstream
from .cloud.base import path_conversion

sentinel = object()
log = logging.getLogger(__name__)


class Fingerprint:
    """Generate a snapshot of the system info."""
    def __init__(self, note=None, id=None, pid=False):
        self.datetime = datetime.utcnow().isoformat()
        self.user = getuser()
        self.version = str(jetstream.__version__)
        self.args = ' '.join(sys.argv)
        self.hostname = gethostname()
        self.pwd = os.getcwd()
        self.note = str(note)
        self.id = id or jetstream.guid()
        if pid:
            self.pid = os.getpid()

    def to_dict(self):
        return vars(self)

    def to_yaml(self):
        return dumps_yaml(vars(self))

    def to_json(self):
        return dumps_json(vars(self), sort_keys=True)


def coerce_tuple(obj):
    """Coerce an object to a tuple.
    Since strings are sequences, iterating over an object that can be a string
    or a sequence can be difficult. This solves the problem by ensuring scalars
    are converted to sequences. """
    if obj is None:
        return tuple()
    elif isinstance(obj, str):
        return (obj, )
    else:
        return tuple(obj)


def coerce_list(obj):
    """Coerce an object to a list.
    Since strings are sequences, iterating over an object that can be a string
    or a sequence can be difficult. This solves the problem by ensuring scalars
    are converted to sequences. """
    if obj is None:
        return list()
    elif isinstance(obj, str):
        return [obj, ]
    else:
        return list(obj)


def config_stack(*sources):
    """Uses Confuse config flattener to merge together nested config sources.
    The main goal here is to allow multiple sources to merge smoothly rather
    than completely overwrite contents for nested properties.

    In this example, a is the project config and b is the config from command
    arguments. We want the command arguments to overwrite the project config,
    but a standard a.update(b) would not give the desired results

        a = {'foo': {'bar': 42, 'baz': 24}}
        b = {'foo': {'bar': 123}}

        >>> a.update(b)
        >>> a
        {'foo': {'bar': 123}}

        Instead config stack gives us a union:

        >>> config_stack(a, b)
        {'foo': {'bar': 123, 'baz': 24}}


    Warning: this is accomplished with a couple hacks -

    1) To prevent Confuse from searching the filesystem for existing config
    files, an unlikely name is added to the Configuration object when it's
    created.

    2) The resulting config object is json dumped then loaded to get convert
    the ordered dictionaries to standard dicts. So.. high performance is not a
    focus here.

    3) Appending/extending a nested array is still not possible:

        >>> b = {'foo': {'bar': [42, 24,]}}
        >>> a = {'foo': {'bar': [123,]}}
        >>> jetstream.utils.config_stack(a, b)
        {'foo': {'bar': [42, 24]}}

    """
    # Hack to prevent Confuse from accidentally finding a config file on the
    # system when were just trying to create an in-memory config object
    n = 'UNLIKELYAPPNAME-' + str(uuid4())
    conf = confuse.Configuration(n, read=False)
    for s in sources:
        if s:
            if isinstance(s, Mapping):
                conf.set(confuse.ConfigSource(s))
            else:
                err = 'Config sources should return Mapping-type objects'
                raise ValueError(err)

    # Hack to remove the ordereddicts, they're just ugly to look at and when
    # working with primitive data types coming from json/yaml files they're
    # generally useless.
    return json.loads(json.dumps(conf.flatten()))


def dict_update_dot_notation(d, key, value):
    """Updates a keys in dictionaries and allows access to nested dictionaries
    with a dot notation:
    Example:
        >>> d = {'foo': {'bar': 'valA', 'baz': 'valb'}}
        >>> dict_update_dot_notation(d, 'foo.baz', 42)
        >>> d
        {'foo': {'bar': 'valA', 'baz': 42}}
    """
    path = key.split('.')
    k = path.pop(0)

    while 1:
        try:
            nk = path.pop(0)
        except IndexError:
            break

        try:
            d = d[k]
        except KeyError:
            d[k] = dict()
            d = d[k]

        k = nk

    d[k] = value


def dict_lookup_dot_notation(d, path):
    path = path.split('.')

    k = path.pop(0)
    v = d[k]
    d = v

    while path:
        k = path.pop(0)
        v = d[k]
        d = v

    return v


def dump_json(obj, *args, sort_keys=True, **kwargs):
    """Attempt to dump `obj` to a JSON file"""
    return json.dump(obj, sort_keys=sort_keys, *args, **kwargs)


def dumps_json(obj, *args, sort_keys=True, **kwargs):
    """Attempt to dump `obj` to a JSON string"""
    return json.dumps(obj, sort_keys=sort_keys, *args, **kwargs)


def dump_yaml(obj, stream):
    """Attempt to dump `obj` to a YAML file"""
    return yaml.dump(obj, stream=stream, default_flow_style=False)


def dumps_yaml(obj):
    """Attempt to dump `obj` to a YAML string"""
    stream = io.StringIO()
    yaml.dump(obj, stream=stream, default_flow_style=False)
    return stream.getvalue()


def dynamic_import(path):
    """Imports functions from other libraries when needed"""
    m, _, f = path.rpartition('.')

    try:
        if m and jetstream.settings['dynamic_import'].get():
            mod = importlib.import_module(m)
            return getattr(mod, f)
        else:
            return getattr(builtins, f)
    except AttributeError as e:
        err = f'Failed to import "{path}": Error: {e}'
        raise AttributeError(err) from None


def guess_local_cpus( default=1):
    return cpu_count() or default


def guess_max_forks(default=500):
    """Returns 1/4 of current ulimit -u value. This leaves a good amount of
    room for subprocesses to continue."""
    try:
        res = int(0.25 * int(subprocess.check_output('ulimit -u', shell=True)))
        return res
    except subprocess.CalledProcessError as e:
        log.debug('Guessing max forks with ulimit -u failed, using default')
        return default


def is_gzip(path, magic_number=b'\x1f\x8b'):
    """Returns True if the path is gzipped."""
    if os.path.exists(path) and not os.path.isfile(path):
        err = 'This should only be used with regular files because otherwise ' \
              'it will lose some data.'
        raise ValueError(err)

    with open(path, 'rb') as fp:
        if fp.read(2) == magic_number:
            return True
        else:
            return False


def is_multiline(s):
    """Returns True if the string contains more than one line.
    This tends to perform much better than len(s.splitlines()) or s.count('\n')
    :param s: string to check
    :return: Bool
    """
    start = 0
    lines = 1

    while True:
        idx = s.find('\n', start)

        if idx == -1:
            return lines > 1
        else:
            lines += 1

        if lines > 1:
            return True
        else:
            start = idx + 1


def is_scalar(obj):
    """Returns `True` if the `obj` should be considered a scalar for YAML
    serialization. Strings are an instance of Sequence, so this function might
    look odd, but it does the job."""
    if isinstance(obj, (str,)):
        return True

    if isinstance(obj, (Sequence, Mapping)):
        return False
    else:
        return True


def filter_records(records, criteria):
    """Given a list of mapping objects (`records`) and a criteria mapping,
    this function returns a list of objects that match filter criteria."""
    matches = list()
    for i in records:
        for k, v in criteria.items():
            if k not in i:
                log.debug('Drop "{}" due to "{}" not in doc'.format(i, k))
                break
            if i[k] != v:
                log.debug('Drop "{}" due to "{}" not == "{}"'.format(i, k, v))
                break
        else:
            matches.append(i)
    return matches


def find(path, name=None):
    """Similar to BSD find, finds files matching name pattern."""
    for dirname, subdirs, files in os.walk(path):
        for f in files:
            if name is None:
                yield os.path.join(dirname, f)
            elif fnmatch.fnmatch(f, name):
                yield os.path.join(dirname, f)


def parse_bool(data):
    """Parse a string value to bool"""
    if data.lower() in ('yes', 'true',):
        return True
    elif data.lower() in ('no', 'false',):
        return False
    else:
        err = f'"{data}" could not be interpreted as a boolean'
        raise TypeError(err)


def parse_csv(data):
    """Parse csv with headers, returns list of dicts"""
    return parse_table(data, dialect='unix', headers=True)


def parse_csv_nh(data):
    """Parse csv with no header, returns list of lists"""
    return parse_table(data, dialect='unix', headers=False)


def parse_json(data):
    return json.loads(data)


def parse_table(data, dialect='unix', headers=True, ordered=False):
    """Attempts to load a table file in any format. Returns a list of 
    dictionaries (or list of lists if no header is available). This requires 
    does not handle comment lines."""    
    if headers:
        rows = csv.DictReader(data.splitlines(), dialect=dialect)
        if not ordered:
            rows = [dict(r) for r in rows]
    else:
        rows = csv.reader(data.splitlines(), dialect=dialect)

    return list(rows)

    
def parse_tsv(data):
    """Parse tsv with headers, returns list of dicts"""
    return parse_table(data, dialect='excel-tab', headers=True)


def parse_tsv_nh(data):
    """Parse tsv with no header, returns list of lists"""
    return parse_table(data, dialect='excel-tab', headers=False)


def parse_txt(data):
    """Parse txt data, returns list of lines"""
    return data.splitlines()


def parse_yaml(data):
    import math
    # import yaml.scanner.ScannerError
    # zfill_val = math.ceil(math.log(len(data.split('\n')), 10))
    # for line_no, line in enumerate(data.split('\n'), start=1):
    #     print('{}  {}'.format(str(line_no).zfill(zfill_val), line))
    try:
        return yaml.safe_load(data)
    except Exception as e:
        except_line = e.problem_mark.line + 1
        log.info('Except line in {}'.format(except_line))
        show_range = set(range(except_line - 50, except_line + 25))
        for line_no, line in enumerate(data.split('\n'), start=1):
            if line_no in show_range:
                print('{}{}  {}'.format('>' if line_no == except_line else ' ', line_no, line))
        raise


def _load(path):
    # todo, if path is valid uri, download
    with open(path, 'r') as fp:
        return fp.read()


def load_file(path, filetype=None):
    """Attempts to load a data file from path, raises :ValueError
    if an suitable loader function is not found in loaders"""
    if filetype is None:
        for ext, fn in loaders.items():
            if path.endswith(ext):
                return fn(path)
        else:
            err = f'No loader extension pattern matched path: {path}'
            raise ValueError(err)
    else:
        try:
            return loaders[filetype](path)
        except KeyError:
            err = f'No loader defined for {filetype}'
            raise ValueError(err)


def load_tsv(path):
    """Load tsv with headers, returns list of dicts"""
    return parse_tsv(_load(path))


def load_csv(path):
    """Load csv with headers, returns list of dicts"""
    return parse_csv(_load(path))


def load_tsv_nh(path):
    """Load tsv with no header, returns list of lists"""
    return parse_tsv_nh(_load(path))


def load_csv_nh(path):
    """Load csv with no header, returns list of lists"""
    return parse_csv_nh(_load(path))


def load_json(path):
    """Load a json file from path"""
    return parse_json(_load(path))


def load_txt(path):
    """Load plain text file as list of lines"""
    return parse_txt(_load(path))


def load_yaml(path):
    """Load a yaml file from `path`"""
    return parse_yaml(_load(path))


def read_lines_allow_gzip(path):
    """Reads line-separated text files, handles gzipped files and recognizes
    universal newlines. This can cause bytes to be lost when reading from a
    pipe. """
    if is_gzip(path):
        with gzip.open(path, 'rb') as fp:
            data = fp.read().decode('utf-8')
    else:
        with open(path, 'r') as fp:
            data = fp.read()
    lines = data.splitlines()
    return lines


def records_to_csv(records, outpath):
    """Writes records (list of dictionaries) out to a csv file"""
    if os.path.exists(outpath):
        raise FileExistsError(outpath)

    keys = set()
    for record in records:
        keys = keys.union(set(record.keys()))

    log.debug('Found keys: {}'.format(keys))

    with open(outpath, 'w') as fp:
        dw = csv.DictWriter(fp, keys)
        dw.writeheader()
        for record in records:
            # Convert non-scalar values to json strings
            for key, value in record.items():
                if not is_scalar(value):
                    record[key] = json.dumps(value)
        dw.writerows(records)


def remove_prefix(string, prefix):
    if string.startswith(prefix):
        return string[len(prefix):]
    else:
        return string


loaders = {k: dynamic_import(v) for k, v in jetstream.settings['loaders'].get(dict).items()}
parsers = {k: dynamic_import(v) for k, v in jetstream.settings['parsers'].get(dict).items()}

# Cloud utils
def replace_root_dir(path, root_dir):
    if not root_dir.endswith(os.path.sep):
        root_dir = root_dir + os.path.sep
    if f'.{os.path.sep}' in path:
        return root_dir + path.split(f'.{os.path.sep}')[1]
    else:
        if path.startswith(os.path.sep):
            return path
        else:
            return root_dir + os.path.basename(path) # path

# Inputs: Local path and remote working directotry
# Outputs: Local to remote mapping
def map_stage_in(lpath, rwd):
    if '->' in lpath:
        return lpath

    rpath = replace_root_dir(lpath, rwd)
    return '{} -> {}'.format(lpath, rpath)

# Inputs: Remote path and local directory
# Outputs: Local to remote mapping
def map_stage_out(rpath, lwd = None):
    if '<-' in rpath:
        return rpath

    lpath = replace_root_dir(rpath, lwd or os.getcwd())
    os.makedirs(os.path.dirname(lpath), exist_ok=True)
    return '{} <- {}'.format(lpath, rpath)

def construct_cjs_cmd(task_body, service_url, cloud_storage=None, cjs_stagein=None, cjs_stageout=None, cloud_downloads=None,
                      cloud_uploads=None, stdout=None, stderr=None, redirected=False, rwd=None, cloud_scripts_dir='.',
                      task_name='', singularity_container_uri=None):
    """
    cloud_downloads will be a blob path
    cloud_uploads will be remote paths, but put on cloud storage with the full relative path as the name
    """
    rwd = rwd or '/tmp/pworks'
    cjs_stagein = cjs_stagein or list()
    cjs_stageout = cjs_stageout or list()
    
    # Create scripts to download/upload data to/from the remote worker
    if cloud_storage is not None:
        # Create script to download data inputs onto the remote node
        remote_download_cmd = cloud_storage.remote_download_cmd(cloud_downloads)
        remote_download_sh_path = os.path.join(cloud_scripts_dir, '.', f'{task_name}.cloud_download.sh')
        with open(remote_download_sh_path, 'w') as down_out:
            down_out.write(remote_download_cmd)
        cjs_stagein.append(remote_download_sh_path)
        
        # Create script to upload data outputs from the remote node into cloud storage
        remote_upload_cmd = cloud_storage.remote_upload_cmd(cloud_uploads)
        remote_upload_sh_path = os.path.join(cloud_scripts_dir, '.', f'{task_name}.cloud_upload.sh')
        with open(remote_upload_sh_path, 'w') as up_out:
            up_out.write(remote_upload_cmd)
        cjs_stagein.append(remote_upload_sh_path)
    
    # Create script to run the main task body
    log.debug('Writing cloud script for {}'.format(task_name))
    cmd_sh_path = os.path.join(cloud_scripts_dir, './{}.sh'.format(task_name))
    with open(cmd_sh_path, 'w') as cmd_sh_out:
        cmd_sh_out.write(task_body)
    cjs_stagein.append(cmd_sh_path)
    
    # Append hostname out
    hostname_out_path = f'./.{task_name}.hostname'
    cjs_stageout.append(hostname_out_path)
    
    # Set paths for remote stdout/stderr
    remote_stdout_path = f'./{task_name}.remote.out'
    remote_stderr_path = f'./{task_name}.remote.err'
    cjs_stageout.append(remote_stdout_path)
    cjs_stageout.append(remote_stderr_path)
    
    # Fill template for stagein and stageout arguments
    input_maps = ' : '.join([map_stage_in(inp, rwd) for inp in cjs_stagein])
    if input_maps:
        input_maps = ' -stagein "{}"'.format(input_maps)
    output_maps = ' : '.join([map_stage_out(outp) for outp in cjs_stageout])
    if output_maps:
        output_maps = ' -stageout "{}"'.format(output_maps)
    
    # Fill in template for redirected, stdout, stderr arguments
    std = '{redirected} {stdout} {stderr}'.format(
        redirected='-redirected' if redirected else '',
        stdout=' -stdout "{}"'.format(stdout) if stdout is not None else '',
        stderr=' -stderr "{}"'.format(stderr) if stderr is not None else ''
    )
    
    # Fill in template to execute bash scripts on the worker node
    if singularity_container_uri is None:
        task_body_cmd = 'bash {}.sh >>{} 2>>{}'.format(task_name, remote_stdout_path, remote_stderr_path) 
    else:
        task_body_cmd = 'singularity exec --cleanenv --nv {container_uri} bash {task_name}.sh >>{remote_stdout_path} 2>>{remote_stderr_path}'.format(
            container_uri=path_conversion(singularity_container_uri),
            task_name=task_name,
            remote_stdout_path=remote_stdout_path,
            remote_stderr_path=remote_stderr_path
        )

    if cloud_storage is not None:
        remote_download_cmd = 'bash {} >>{} 2>>{};'.format(os.path.basename(remote_download_sh_path), remote_stdout_path, remote_stderr_path)
        remote_upload_cmd = 'bash {} >>{} 2>>{};'.format(os.path.basename(remote_upload_sh_path), remote_stdout_path, remote_stderr_path)
    else:
        remote_download_cmd, remote_upload_cmd = '', ''
    
    # Fill in template for complete cog-job-submit command
    cjs_final_cmd = (
        'cog-job-submit -provider "coaster-persistent" -attributes "maxWallTime=240:00:00" '
        '{std} -service-contact "{service_url}"{input_maps}{output_maps} -directory "{rwd}" '
        '/bin/bash -c "mkdir -p {rwd};cd {rwd};hostname > {hostname_log};'
        '{cloud_download_cmd}{cmd};{cloud_upload_cmd}"'
    ).format(
        std=std,
        service_url=service_url,
        input_maps=input_maps,
        output_maps=output_maps,
        rwd=rwd,
        cmd=task_body_cmd,
        cloud_download_cmd=remote_download_cmd,
        cloud_upload_cmd=remote_upload_cmd,
        hostname_log=hostname_out_path
    )
    
    try:
        with open(os.path.join(cloud_scripts_dir, f'{task_name}.cjs'), 'w') as cjs_final_out:
            cjs_final_out.write(cjs_final_cmd + '\n')
    except:
        pass  # Fail silently
    
    return cjs_final_cmd