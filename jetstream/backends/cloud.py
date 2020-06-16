import logging
import asyncio
import jetstream
from jetstream.backends import BaseBackend
from ..cloud.azure import AzureStorageSession
from asyncio import BoundedSemaphore, create_subprocess_shell, CancelledError
import os
import sys
import requests
import json
from random import randint
import re
from time import sleep
import sys
import tempfile
import glob
import subprocess
import re
from datetime import datetime
import urllib
import uuid

log = logging.getLogger('jetstream.local')
# cjs_dir_path = "/home/avidalto/projects/2020/jetstream/jetstream-may/cjs-backend/"
# api_key = "649cec38be4772725c5c06d4b0e4d494"
# pool_name = "linuxpool"

def get_pool_info(pool_name, api_key, pworks_url="http://beta.parallel.works"):
    rq = requests.get('{pworks_url}/api/resources?key={api_key}'.format(pworks_url=pworks_url, api_key=api_key))
    for pool_data in rq.json():
        if pool_data['name'] == pool_name:
            log.info('Got pooldata: {}'.format(pool_data))
            return {
                'serviceport': pool_data['info']['ports']['serviceport'],
                'controlport': pool_data['info']['ports']['controlport'],
                'maxworkers': int(pool_data['settings']['max']),
                # 'cpus': int(pool_data['info']['cpuPerWorker']) // 2   This is temporary
                'cpus': 8
            }


# Same function as in workflow.py jetstream script
def compile(pattern):
    return re.compile('^{}$'.format(pattern))

# Finds the files that match a regex pattern
def expand_regex_path(patterns):
    pfiles = []
    for pattern in patterns:
        if pattern.startswith("~"):
            pattern = os.path.expanduser("~") + pattern[1:]
        if pattern.startswith("/"):
            rootdir = ""
            for subdir in pattern.split("/")[1:]:
                if os.path.isdir(rootdir + "/" + subdir):
                    rootdir = rootdir + "/" + subdir
                else:
                    break
        else:
            rootdir = "./"
        cpattern = compile(pattern)
        for root, dirs, files in os.walk(rootdir):
            for f in files:
                path = root + "/" + f
                if rootdir == "./" and pattern[0:2] != "./":
                    path = path[2:]
                if cpattern.match(path):
                    pfiles.append(path)
    return pfiles


def parse_reference_input(ref_input_directive):
    ref_input = list()
    if isinstance(ref_input_directive, str):
        ref_input_directive = [ref_input_directive]
    for ref in ref_input_directive:
        ref_input.extend(glob.glob(ref))
    return ref_input


def get_cloud_directive(key, task_directives, cloud_args_key='cloud_args'):
    return task_directives.get(cloud_args_key, dict()).get(key)


class CloudSwiftBackend(BaseBackend):
    def __init__(self, pool_name=None, api_key=None, cpus=None, blocking_io_penalty=None, wrapper_frontmatter=None,
                 **kwargs):
        """The LocalBackend executes tasks as processes on the local machine.

        This contains a semaphore that limits tasks by the number of cpus
        that they require. It requires that self.runner be set to get the
        event loop, so it's not instantiated until preflight.

        :param cpus: If this is None, the number of available CPUs will be
            guessed. This cannot be changed after starting the backend.
        :param blocking_io_penalty: Delay (in seconds) when a BlockingIOError
            prevents a new process from spawning.
        :param max_concurrency: Max concurrency limit
        """
        super().__init__()
        self.pool_info = get_pool_info(pool_name, api_key)
        log.info('Pool info: {}'.format(self.pool_info))
        self.cpus = self.pool_info['cpus'] * self.pool_info['maxworkers']
        self.bip = blocking_io_penalty \
                   or jetstream.settings['backends']['local']['blocking_io_penalty'].get(int)
        self._cpu_sem = BoundedSemaphore(int(self.cpus))
        self.project_dir = os.getcwd()
        try:
            os.remove('cjs_cmds_debug.sh')
        except:
            pass  # Fail silently
        
        # Make directory for cjs launch scripts
        self.cloud_scripts_dir = os.path.join(self.project_dir, 'cloud_scripts')
        os.makedirs(self.cloud_scripts_dir, exist_ok=True)
        
        # Store default config wrapper frontmatter
        self.wrapper_frontmatter = wrapper_frontmatter
        
        # If we need it, store a temporary name for a cloud container
        self.cloud_storage = AzureStorageSession(
            **kwargs['azure_params'],
            create_temp_container=False
        )
        self.cloud_storage._temp_container_name = 'test'
        
        log.info(f'CloudSwiftBackend initialized with {self.cpus} cpus')

    async def spawn(self, task):
        is_local_task = task.directives.get('cloud_args', dict()).get('local_task', False)
        log.info('Spawn ({}): {}'.format('Local' if is_local_task else 'Cloud', task))
        
        if is_local_task:
            return await self.spawn_local(task)
        
        # This is a cloud task
        return await self.spawn_cloud(task)

    async def spawn_cloud(self, task):
        if 'cmd' not in task.directives:
            return task.complete()
        
        cmd = task.directives['cmd'] 
        cpus = task.directives.get('cpus', 0)
        cpus_reserved = 0
        open_fps = list()

        if cpus > self.pool_info['cpus']:
            #raise RuntimeError('Task cpus greater than available cpus')
            raise RuntimeError(
                'Task cpus ({}) greater than available cpus ({}) in worker'.format(cpus, self.pool_info['cpus'])
            )
        
        try:
            # Get file descriptor paths and file pointers for this task
            fd_paths = {
                fd_name: fd
                for fd_name, fd in zip(('stdin', 'stdout', 'stderr'), self.get_fd_paths(task))
            }
            fd_filepointers = {
                fd_name: open(fd, fd_mode) if fd else None
                for fd_mode, (fd_name, fd) in zip(('r', 'w', 'w'), fd_paths.items())
            }
            
            # Upload data inputs into cloud storage
            log.info('Making call to input data into the cloud')
            blob_inputs_to_remote(task.directives['input'], self.cloud_storage)
            
            # Upload reference inputs into cloud storage
            reference_inputs = parse_reference_input(task.directives.get('cloud_args', dict()).get('reference_input', list()))
            blob_inputs_to_remote(reference_inputs, self.cloud_storage, blob_basename=True)        
            
            # If the user provides a non-URL path to a container and explicitly says it should be transfer,
            # then consider it similar to reference data and upload it to cloud storage
            singularity_container_uri = task.directives.get('cloud_args', dict()).get('singularity_container')
            container_input = (
                [singularity_container_uri]
                if (
                    singularity_container_uri is not None
                    and not urllib.parse.urlparse(singularity_container_uri).scheme
                    and get_cloud_directive('transfer_container_to_remote', task.directives)
                )
                else list()
            )
            blob_inputs_to_remote(container_input, self.cloud_storage)

            # Construct the cog-job-submit command for execution
            cjs_cmd = construct_cjs_cmd(
                task_body=cmd,
                service_url='http://beta.parallel.works:{}'.format(self.pool_info['serviceport']),
                cjs_stagein=None,
                cjs_stageout=None,
                cloud_downloads=task.directives['input'] + reference_inputs + container_input,
                cloud_uploads=task.directives['output'],
                account_name=self.cloud_storage.storage_account_name,
                account_key=self.cloud_storage.storage_account_key,
                cloud_scripts_dir=self.cloud_scripts_dir,
                container=self.cloud_storage._temp_container_name,
                singularity_container_uri=singularity_container_uri,
                task_name=task.name
            )
            
            # Async submit as a subprocess
            p = await self.subprocess_sh(
                cjs_cmd,
                **fd_filepointers
            )

            # Once command is executed, update task with some process metadata
            task.state.update(
                stdout_path=fd_paths['stdout'],
                stderr_path=fd_paths['stderr'],
                label=f'CloudSwift({p.pid})',
            )

            log.info(f'CloudSwiftBackend spawned({p.pid}): {task.name}')
            rc = await p.wait()

            if rc != 0:
                log.info(f'Failed: {task.name}')
                return task.fail(p.returncode)
            else:
                # Download completed data files from cloud storage
                blob_outputs_to_local(task.directives['output'], self.cloud_storage)
                log.info(f'Complete: {task.name}')
                return task.complete(p.returncode)
        except CancelledError:
            task.state['err'] = 'Runner cancelled Backend.spawn()'
            return task.fail(-15)
        except Exception as e:
            log.error('Exception: {}'.format(e))
            import traceback
            traceback.print_exc()
            raise
        finally:
            for fp in fd_filepointers.values():
                if fp is not None:
                    fp.close()

            for i in range(cpus_reserved):
                self._cpu_sem.release()

            return task
    
    async def spawn_local(self, task):
        if 'cmd' not in task.directives:
            return task.complete()

        cmd = task.directives['cmd']
        cpus = task.directives.get('cpus', 0)
        cpus_reserved = 0
        open_fps = list()

        if cpus > self.cpus:
            raise RuntimeError('Task cpus greater than available cpus')

        try:
            for i in range(task.directives.get('cpus', 0)):
                await self._cpu_sem.acquire()
                cpus_reserved += 1

            stdin, stdout, stderr = self.get_fd_paths(task)

            if stdin:
                stdin_fp = open(stdin, 'r')
                open_fps.append(stdin_fp)
            else:
                stdin_fp = None

            if stdout:
                stdout_fp = open(stdout, 'w')
                open_fps.append(stdout_fp)
            else:
                stdout_fp = None

            if stderr:
                stderr_fp = open(stderr, 'w')
                open_fps.append(stderr_fp)
            else:
                stderr_fp = None

            p = await self.subprocess_sh(
                cmd,
                stdin=stdin_fp,
                stdout=stdout_fp,
                stderr=stderr_fp
            )

            task.state.update(
                stdout_path=stdout,
                stderr_path=stderr,
                label=f'Slurm({p.pid})',
            )

            log.info(f'LocalBackend spawned({p.pid}): {task.name}')
            rc = await p.wait()

            if rc != 0:
                log.info(f'Failed: {task.name}')
                return task.fail(p.returncode)
            else:
                log.info(f'Complete: {task.name}')
                return task.complete(p.returncode)
        except CancelledError:
            task.state['err'] = 'Runner cancelled Backend.spawn()'
            return task.fail(-15)
        finally:
            for fp in open_fps:
                fp.close()

            for i in range(cpus_reserved):
                self._cpu_sem.release()

            return task

    async def subprocess_sh(
            self, args, *, stdin=None, stdout=None, stderr=None,
            cwd=None, encoding=None, errors=None, env=None,
            loop=None, executable="/bin/bash"):
        """Asynchronous version of subprocess.run

        This will always use a shell to launch the subprocess, and it prefers
        /bin/bash (can be changed via arguments)"""
        log.debug(f'subprocess_sh:\n{args}')
        #print("{} {}".format(executable, args))
        while 1:
            try:
                p = await create_subprocess_shell(
                    args,
                    stdin=stdin,
                    stdout=stdout,
                    stderr=stderr,
                    cwd=cwd,
                    encoding=encoding,
                    errors=errors,
                    env=env,
                    loop=loop,
                    executable=executable
                )
                break
            except BlockingIOError as e:
                log.warning(f'System refusing new processes: {e}')
                await asyncio.sleep(self.bip)

        return p


def replace_root_dir(path, root_dir):
    if not root_dir.endswith("/"):
        root_dir = root_dir + "/"
    if "./" in path:
        return root_dir + path.split("./")[1]
    else:
        if path.startswith("/"):
            return path
        else:
            return root_dir + os.path.basename(path) # path

# Inputs: Local path and remote working directotry
# Outputs: Local to remote mapping
def map_stage_in(lpath, rwd):
    if "->" in lpath:
        return lpath

    rpath = replace_root_dir(lpath, rwd)
    return "{} -> {}".format(lpath, rpath)

# Inputs: Remote path and local directory
# Outputs: Local to remote mapping
def map_stage_out(rpath, lwd = None):
    if "<-" in rpath:
        return rpath

    if lwd is None:
        lwd = os.getcwd()
    lpath = replace_root_dir(rpath, lwd)
    os.makedirs(os.path.dirname(lpath), exist_ok=True)
    return "{} <- {}".format(lpath, rpath)


def dummy_dot_in_path(path):
    return '.' in path.split(os.path.sep)
        

def path_conversion(path):
    if dummy_dot_in_path(path):
        # local: /absolute/path/./to/file ----> remote: current/dir/to/file
        # local: relative/path/./to/file ----> remote: current/dir/to/file
        # local: current/dir/to/file <---- remote: /absolute/path/./to/file
        # local: current/dir/to/file <---- remote: relative/path/./to/file
        return path.rsplit(f'.{os.path.sep}')[-1]
    
    if os.path.isabs(path):
        # local: /absolute/path/to/file ----> remote: /absolute/path/to/file
        # local: /absolute/path/to/file <---- remote: /absolute/path/to/file
        return path
    
    # local: relative/path/to/file ----> remote: current/dir/file
    # local: current/dir/file <---- remote: relative/path/to/file
    return path.split(os.path.sep)[-1]
        

def blob_inputs_to_remote(inputs, cloud_storage, container=None, blob_basename=False):
    for input_file in inputs:
        cloud_storage.upload_blob(
            input_file,
            blobpath=os.path.basename(input_file) if blob_basename else path_conversion(input_file)
        )
    

def blob_outputs_to_local(outputs, cloud_storage, container=None, blob_basename=False):
    for output_file in outputs:
        try:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            cloud_storage.download_blob(
                output_file,
                blobpath=os.path.basename(output_file) if blob_basename else path_conversion(output_file)
            )
        except Exception as e:
            log.info(f'Error in download: {e}')
    

def construct_cjs_cmd(task_body, service_url, cjs_stagein=None, cjs_stageout=None, cloud_downloads=None, cloud_uploads=None,
                       stdout=None, stderr=None, redirected=False, rwd=None, container='temp', account_name='', account_key='',
                       cloud_scripts_dir='', task_name='', singularity_container_uri=None):
    """
    cloud_downloads will be a blob path
    cloud_uploads will be remote paths, but put on cloud storage with the full relative path as the name
    """
    rwd = rwd or '/tmp/pworks/{}'.format(str(randint(0, 99999)).zfill(5))
    
    # Create script to download data inputs onto the remote node
    cloud_download_cmd_template = (
        'if [[ ! -f "{blobpath}" ]];then mkdir -p {blobpath_dirname}; '
        'singularity exec {az_sif_path} az storage blob download --name {blobpath} --file {blobpath} '
        '--container-name {container_name} --account-name {account_name} --account-key {account_key};fi\n'
    )
    cloud_download_cmd = ''
    for cloud_download in cloud_downloads or list():
        cloud_download = path_conversion(cloud_download)
        cloud_download_cmd += cloud_download_cmd_template.format(
            blobpath=cloud_download,
            blobpath_dirname=os.path.dirname(cloud_download),
            az_sif_path='docker://mcr.microsoft.com/azure-cli',
            container_name=container,
            account_name=account_name,
            account_key=account_key
        )
    cloud_download_sh_path = os.path.join(cloud_scripts_dir, '.', 'cloud_download_{}.sh'.format(str(uuid.uuid4())[:8]))
    with open(cloud_download_sh_path, 'w') as down_out:
        down_out.write(cloud_download_cmd)
    
    # Create script to upload data outputs from the remote node into cloud storage
    cloud_upload_cmd_template = (
        'singularity exec {az_sif_path} az storage blob upload --name {blobpath} --file {blobpath} '
        '--container-name {container_name} --account-name {account_name} --account-key {account_key};\n'
    )
    cloud_upload_cmd = ''
    for cloud_upload in cloud_uploads or list():
        cloud_upload = path_conversion(cloud_upload)
        cloud_upload_cmd += cloud_upload_cmd_template.format(
            blobpath=cloud_upload,
            blobpath_dirname=os.path.dirname(cloud_upload),
            az_sif_path='docker://mcr.microsoft.com/azure-cli',
            container_name=container,
            account_name=account_name,
            account_key=account_key
        )
    cloud_upload_sh_path = os.path.join(cloud_scripts_dir, '.', 'cloud_upload_{}.sh'.format(str(uuid.uuid4())[:8]))
    with open(cloud_upload_sh_path, 'w') as up_out:
        up_out.write(cloud_upload_cmd)
    
    # Create script to run the main task body
    log.debug('Writing cloud script for {}'.format(task_name))
    cmd_sh_path = os.path.join(cloud_scripts_dir, './{}.sh'.format(task_name))
    with open(cmd_sh_path, 'w') as cmd_sh_out:
        cmd_sh_out.write(task_body)
        
    # Append the three above scripts to stagein list
    cjs_stagein = cjs_stagein or list()
    cjs_stagein.append(cloud_download_sh_path)
    cjs_stagein.append(cloud_upload_sh_path)
    cjs_stagein.append(cmd_sh_path)
    
    # Fill template for stagein and stageout arguments
    input_maps = ' : '.join([map_stage_in(inp, rwd) for inp in (cjs_stagein or list())])
    if input_maps:
        input_maps = ' -stagein "{}"'.format(input_maps)
    output_maps = ' : '.join([map_stage_out(outp) for outp in (cjs_stageout or list())])
    if output_maps:
        output_maps = ' -stageout "{}"'.format(output_maps)
    
    # Fill in template for redirected, stdout, stderr arguments
    std = '{redirected} {stdout} {stderr}'.format(
        redirected='-redirected' if redirected else '',
        stdout=' -stdout "{}"'.format(stdout) if stdout is not None else '',
        stderr=' -stderr "{}"'.format(stderr) if stderr is not None else ''
    )
    
    # Fill in template to execute bash scripts on the worker node
    task_body_cmd = (
        'bash {}.sh'.format(task_name) if singularity_container_uri is None
        else 'singularity exec --cleanenv --nv {container_uri} bash {task_name}.sh'.format(
            container_uri=singularity_container_uri,
            task_name=task_name
        )
    )
    cloud_download_cmd = 'bash {};'.format(os.path.basename(cloud_download_sh_path))
    cloud_upload_cmd = 'bash {};'.format(os.path.basename(cloud_upload_sh_path))
    
    # Fill in template for complete cog-job-submit command
    return (
        'cog-job-submit -provider "coaster-persistent" -attributes "maxWallTime=240:00:00" '
        '{std} -service-contact "{service_url}"{input_maps}{output_maps} -directory "{rwd}" '
        '/bin/bash -c "mkdir -p {rwd};cd {rwd};{cloud_download_cmd}{cmd};{cloud_upload_cmd}"'
    ).format(
        std=std,
        service_url=service_url,
        input_maps=input_maps,
        output_maps=output_maps,
        rwd=rwd,
        cmd=task_body_cmd,
        cloud_download_cmd=cloud_download_cmd,
        cloud_upload_cmd=cloud_upload_cmd
    )
