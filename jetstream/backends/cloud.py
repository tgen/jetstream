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
import yaml
import time
import traceback

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
                # 'cpus': int(pool_data['info']['cpuPerWorker']) // pool_data['settings']['jobsPerNode']  # This is temporary
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
        if is_remote_uri(ref):
            ref_input.append(ref)
        else:
            ref_input.extend(glob.glob(ref))
    return ref_input


def get_cloud_directive(key, task_directives, cloud_args_key='cloud_args'):
    return task_directives.get(cloud_args_key, dict()).get(key)


class CloudSwiftBackend(BaseBackend):
    def __init__(self, pw_pool_name=None, pw_api_key=None, cpus=None, blocking_io_penalty=None, wrapper_frontmatter=None,
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
        self.is_pw_pool = pw_pool_name is not None
        if pw_pool_name is not None and pw_api_key is not None:
            self.pool_info = get_pool_info(pw_pool_name, pw_api_key)
            log.info('PW Pool info: {}'.format(self.pool_info))
        else:
            self.pool_info = {
                'cpus': kwargs['pool_info']['cpus_per_worker'],
                'maxworkers': kwargs['pool_info']['workers'],
                'serviceurl': kwargs['pool_info']['serviceurl'],
            }
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
        
        self.cloud_logs_dir = os.path.join(self.project_dir, 'cloud_logs')
        os.makedirs(self.cloud_logs_dir, exist_ok=True)
        
        # Store default config wrapper frontmatter
        self.wrapper_frontmatter = wrapper_frontmatter
        
        # If we need it, store a temporary name for a cloud container
        self.cloud_storage = AzureStorageSession(
            **kwargs['azure_params'],
            create_temp_container=True
            # create_temp_container=False
        )
        if kwargs.get('azure_blob_container') is not None:
            self.cloud_storage._temp_container_name = kwargs['azure_blob_container']
        
        # Initialize cloud metrics log
        CloudMetricsLogger.init(at=os.path.join(self.project_dir, 'cloud_metrics_{}.yaml'.format(datetime.now().strftime('%Y%m%d%H%M%S'))))
        
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
        
        start_time = datetime.now()
        bytes_sent_bundle, bytes_received_bundle = list(), list()
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
            data_metrics = blob_inputs_to_remote(task.directives['input'], self.cloud_storage)
            
            # Upload reference inputs into cloud storage
            reference_inputs = parse_reference_input(task.directives.get('cloud_args', dict()).get('reference_input', list()))
            ref_metrics = blob_inputs_to_remote(reference_inputs, self.cloud_storage, blob_basename=True)        
            
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
            container_metrics = blob_inputs_to_remote(container_input, self.cloud_storage)
            
            # Log metrics for input data
            total_input_metrics = data_metrics + ref_metrics + container_metrics
            for m in total_input_metrics:
                bytes_sent_bundle.append(m)
            bytes_sent_bundle.append({
                'name': 'total',
                'size': sum([max(0, t['size']) for t in total_input_metrics]),
                'time': sum([max(0, t['time']) for t in total_input_metrics])
            })

            # Construct the cog-job-submit command for execution
            cjs_cmd = construct_cjs_cmd(
                task_body=cmd,
                service_url='http://beta.parallel.works:{}'.format(self.pool_info['serviceport']) if self.is_pw_pool else self.pool_info['serviceurl'],
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
            log.debug(cjs_cmd)
            
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
                output_metrics = blob_outputs_to_local(task.directives['output'], self.cloud_storage)
                log.info(f'Complete: {task.name}')
                
                # Log metrics for output data
                for m in output_metrics:
                    bytes_received_bundle.append(m)
                bytes_received_bundle.append({
                    'name': 'total',
                    'size': sum([max(0, t['size']) for t in output_metrics]),
                    'time': sum([max(0, t['time']) for t in output_metrics])
                })
                
                return task.complete(p.returncode)
        except CancelledError:
            task.state['err'] = 'Runner cancelled Backend.spawn()'
            return task.fail(-15)
        except Exception as e:
            log.error('Exception: {}'.format(e))
            traceback.print_exc()
            raise
        finally:
            for fp in fd_filepointers.values():
                if fp is not None:
                    fp.close()

            for i in range(cpus_reserved):
                self._cpu_sem.release()
            
            # Get task runtime and which node it ran on
            elapsed_time = datetime.now() - start_time
            try:
                with open(f'.{task.name}.hostname', 'r') as hostname_log:
                    hostname = hostname_log.read().strip()
                subprocess.call(['mv'] +  glob.glob('*.remote.out') + [self.cloud_logs_dir])
                subprocess.call(['mv'] +  glob.glob('*.remote.err') + [self.cloud_logs_dir])
                os.remove(f'.{task.name}.hostname')
            except:
                pass  # Fail silently
            
            CloudMetricsLogger.write_record({
                'task': task.name,
                'start_datetime': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'elapsed_time': str(elapsed_time),
                'end_datetime': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'in_files': bytes_sent_bundle,
                'out_files': bytes_received_bundle,
                'node': hostname
            })

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
    elapsed_times = list()
    for input_file in inputs:
        try:
            start = time.time()
            cloud_storage.upload_blob(
                input_file,
                blobpath=os.path.basename(input_file) if blob_basename else path_conversion(input_file)
            )
            elapsed_times.append({
                'name': input_file,
                'size': os.stat(input_file).st_size,
                'time': time.time() - start
            })
        except:
            elapsed_times.append({'name': input_file, 'size': -1, 'time': -1})
    return elapsed_times
    

def blob_outputs_to_local(outputs, cloud_storage, container=None, blob_basename=False):
    elapsed_times = list()
    for output_file in outputs:
        try:
            start = time.time()
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            cloud_storage.download_blob(
                output_file,
                blobpath=os.path.basename(output_file) if blob_basename else path_conversion(output_file)
            )
            elapsed_times.append({
                'name': output_file,
                'size': os.stat(output_file).st_size,
                'time': time.time() - start
            })
        except Exception as e:
            log.info(f'Error in download: {e}')
            elapsed_times.append({'name': output_file, 'size': -1, 'time': -1})
    return elapsed_times


def is_cloud_bucket_uri(uri):
    protocol = urllib.parse.urlparse(uri).scheme
    return protocol != '' and protocol.strip().lower() in {'s3', 'gs', 'az'}

def is_remote_uri(uri):
    return urllib.parse.urlparse(uri).scheme != ''


def construct_cloud_download_script(cloud_downloads=None):
    # Create script to download data inputs onto the remote node
    cloud_download_cmd_template = (
        'if [[ ! -f "{blobpath}" ]];then mkdir -p {blobpath_dirname}; '
        'singularity exec {az_sif_path} az storage blob download --name {blobpath} --file {blobpath} '
        '--container-name {container_name} --account-name {account_name} --account-key {account_key};fi\n'
    )
    cloud_download_cmd = ''
    for cloud_download in cloud_downloads or list():
        if is_remote_uri(cloud_download):
            # TODO There is currently no mechanism to tell wget where to put the file
            cloud_download_cmd += 'wget {}\n'.format(cloud_download)
        else:
            cloud_download = path_conversion(cloud_download)
            cloud_download_cmd += cloud_download_cmd_template.format(
                blobpath=cloud_download,
                blobpath_dirname=os.path.dirname(cloud_download),
                az_sif_path='docker://mcr.microsoft.com/azure-cli',
                container_name=container,
                account_name=account_name,
                account_key=account_key
            )
    cloud_download_sh_path = os.path.join(cloud_scripts_dir, '.', f'{task_name}.cloud_download.sh')
    with open(cloud_download_sh_path, 'w') as down_out:
        down_out.write(cloud_download_cmd)


def construct_customcoaster_cmd(cloud_downloads=None, cloud_uploads=None, container='temp', account_name='', account_key='',
                                cloud_scripts_dir='.', task_name='', task_body='', stagein=None, stageout=None, 
                                singularity_container_uri=None, rwd=None, runtask_py_path=''):
    rwd = rwd or '/tmp/pworks'
    
    # Create script to download data inputs onto the remote node
    cloud_download_cmd_template = (
        'if [[ ! -f "{blobpath}" ]];then mkdir -p {blobpath_dirname}; '
        'singularity exec {az_sif_path} az storage blob download --name {blobpath} --file {blobpath} '
        '--container-name {container_name} --account-name {account_name} --account-key {account_key};fi\n'
    )
    cloud_download_cmd = ''
    for cloud_download in cloud_downloads or list():
        if is_remote_uri(cloud_download):
            # TODO There is currently no mechanism to tell wget where to put the file
            url_filepath = urllib.parse.urlparse(cloud_download).path
            cloud_download_cmd += 'if [[ ! -f "{}" ]];then wget {};fi\n'.format(os.path.basename(url_filepath), cloud_download)
        else:
            cloud_download = path_conversion(cloud_download)
            cloud_download_cmd += cloud_download_cmd_template.format(
                blobpath=cloud_download,
                blobpath_dirname=os.path.dirname(cloud_download),
                az_sif_path='docker://mcr.microsoft.com/azure-cli',
                container_name=container,
                account_name=account_name,
                account_key=account_key
            )
    cloud_download_sh_path = os.path.join(cloud_scripts_dir, '.', f'{task_name}.cloud_download.sh')
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
    cloud_upload_sh_path = os.path.join(cloud_scripts_dir, '.', f'{task_name}.cloud_upload.sh')
    with open(cloud_upload_sh_path, 'w') as up_out:
        up_out.write(cloud_upload_cmd)
    
    # Create script to run the main task body
    log.debug('Writing cloud script for {}'.format(task_name))
    cmd_sh_path = os.path.join(cloud_scripts_dir, './{}.sh'.format(task_name))
    with open(cmd_sh_path, 'w') as cmd_sh_out:
        cmd_sh_out.write(task_body)
    
    # Append the three above scripts to stagein list
    stagein = stagein or list()
    stagein.append(cloud_download_sh_path)
    stagein.append(cloud_upload_sh_path)
    stagein.append(cmd_sh_path)
    
    
    # Set paths for hostname, remote stdout/stderr
    stageout = stageout or list()
    hostname_out_path = f'./.{task_name}.hostname'
    stageout.append(hostname_out_path)
    remote_stdout_path = f'./{task_name}.remote.out'
    remote_stderr_path = f'./{task_name}.remote.err'
    stageout.append(remote_stdout_path)
    stageout.append(remote_stderr_path)
    
    # Fill in template to execute bash scripts on the worker node
    task_body_cmd = (
        'bash {}.sh >>{} 2>>{}'.format(task_name, remote_stdout_path, remote_stderr_path) if singularity_container_uri is None
        else 'singularity exec --cleanenv --nv {container_uri} bash {task_name}.sh >>{remote_stdout_path} 2>>{remote_stderr_path}'.format(
            container_uri=singularity_container_uri if urllib.parse.urlparse(singularity_container_uri).scheme else path_conversion(singularity_container_uri),
            task_name=task_name,
            remote_stdout_path=remote_stdout_path,
            remote_stderr_path=remote_stderr_path
        )
    )
    cloud_download_cmd = 'bash {} >>{} 2>>{}'.format(os.path.basename(cloud_download_sh_path), remote_stdout_path, remote_stderr_path)
    cloud_upload_cmd = 'bash {} >>{} 2>>{}'.format(os.path.basename(cloud_upload_sh_path), remote_stdout_path, remote_stderr_path)
    
    # Fill in template for complete cog-job-submit command
    cjs_final_cmd = (
        'mkdir -p {rwd};cd {rwd};hostname > {hostname_out_path};'
        '{cloud_download_cmd};{cmd};{cloud_upload_cmd}'
    ).format(
        # std=std,
        # service_url=service_url,
        # input_maps=input_maps,
        # output_maps=output_maps,
        rwd=rwd,
        cmd=task_body_cmd,
        cloud_download_cmd=cloud_download_cmd,
        cloud_upload_cmd=cloud_upload_cmd,
        hostname_out_path=hostname_out_path
    )
    
    stagein_arg = ' '.join(stagein)
    stageout_arg = ' '.join(stageout)
    worker_cmd = (f'python {runtask_py_path} --port 6000 --host localhost '
                  f'--cmd \'{cjs_final_cmd}\' --rundir {rwd} '
                  f'--stagein {stagein_arg} '
                  f'--stageout {stageout_arg} '
                  + '--inputfiles {} '.format(' '.join(cloud_downloads))
                  + '--outputfiles {}'.format(' '.join(cloud_uploads)))
                  # f'--stageout {hostname_out_path} > /dev/null')
    
    with open(os.path.join(cloud_scripts_dir, f'{task_name}.worker.cmd'), 'w') as worker_final_out:
        worker_final_out.write(worker_cmd + '\n')
    
    return worker_cmd
    

def construct_cjs_cmd(task_body, service_url, cjs_stagein=None, cjs_stageout=None, cloud_downloads=None, cloud_uploads=None,
                       stdout=None, stderr=None, redirected=False, rwd=None, container='temp', account_name='', account_key='',
                       cloud_scripts_dir='.', task_name='', singularity_container_uri=None):
    """
    cloud_downloads will be a blob path
    cloud_uploads will be remote paths, but put on cloud storage with the full relative path as the name
    """
    # rwd = rwd or '/tmp/pworks/{}'.format(str(randint(0, 99999)).zfill(5))
    rwd = rwd or '/tmp/pworks'
    
    # Create script to download data inputs onto the remote node
    cloud_download_cmd_template = (
        'if [[ ! -f "{blobpath}" ]];then mkdir -p {blobpath_dirname}; '
        'singularity exec {az_sif_path} az storage blob download --name {blobpath} --file {blobpath} '
        '--container-name {container_name} --account-name {account_name} --account-key {account_key};fi\n'
    )
    cloud_download_cmd = ''
    for cloud_download in cloud_downloads or list():
        if is_remote_uri(cloud_download):
            # TODO There is currently no mechanism to tell wget where to put the file
            url_filepath = urllib.parse.urlparse(cloud_download).path
            cloud_download_cmd += 'if [[ ! -f "{}" ]];then wget {};fi\n'.format(os.path.basename(url_filepath), cloud_download)
        else:
            cloud_download = path_conversion(cloud_download)
            cloud_download_cmd += cloud_download_cmd_template.format(
                blobpath=cloud_download,
                blobpath_dirname=os.path.dirname(cloud_download),
                az_sif_path='docker://mcr.microsoft.com/azure-cli',
                container_name=container,
                account_name=account_name,
                account_key=account_key
            )
    cloud_download_sh_path = os.path.join(cloud_scripts_dir, '.', f'{task_name}.cloud_download.sh')
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
    cloud_upload_sh_path = os.path.join(cloud_scripts_dir, '.', f'{task_name}.cloud_upload.sh')
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
    
    # Append hostname out
    hostname_out_path = f'./.{task_name}.hostname'
    cjs_stageout = cjs_stageout or list()
    cjs_stageout.append(hostname_out_path)
    
    # Set paths for remote stdout/stderr
    remote_stdout_path = f'./{task_name}.remote.out'
    remote_stderr_path = f'./{task_name}.remote.err'
    cjs_stageout.append(remote_stdout_path)
    cjs_stageout.append(remote_stderr_path)
    
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
        'bash {}.sh >>{} 2>>{}'.format(task_name, remote_stdout_path, remote_stderr_path) if singularity_container_uri is None
        else 'singularity exec --cleanenv --nv {container_uri} bash {task_name}.sh >>{remote_stdout_path} 2>>{remote_stderr_path}'.format(
            container_uri=singularity_container_uri if urllib.parse.urlparse(singularity_container_uri).scheme else path_conversion(singularity_container_uri),
            task_name=task_name,
            remote_stdout_path=remote_stdout_path,
            remote_stderr_path=remote_stderr_path
        )
    )
    cloud_download_cmd = 'bash {} >>{} 2>>{};'.format(os.path.basename(cloud_download_sh_path), remote_stdout_path, remote_stderr_path)
    cloud_upload_cmd = 'bash {} >>{} 2>>{};'.format(os.path.basename(cloud_upload_sh_path), remote_stdout_path, remote_stderr_path)
    
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
        cloud_download_cmd=cloud_download_cmd,
        cloud_upload_cmd=cloud_upload_cmd,
        hostname_log=hostname_out_path
    )
    
    with open(os.path.join(cloud_scripts_dir, f'{task_name}.cjs'), 'w') as cjs_final_out:
        cjs_final_out.write(cjs_final_cmd + '\n')
    
    return cjs_final_cmd


class CloudMetricsLogger:
    _metrics_file = None
    
    @classmethod
    def init(cls, at=None):
        if cls._metrics_file is None:
            _metrics_path = at or 'cloud_metrics_{}.yaml'.format(datetime.now().strftime('%Y%m%d%H%M%S'))
            with open(_metrics_path, 'w') as cloud_metrics:
                cloud_metrics.write(yaml.dump({
                    'run_started': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    'tasks': list()
                }))
            cls._metrics_file = _metrics_path
        
    
    @classmethod
    def write_record(cls, record):
        if cls._metrics_file is not None:
            with open(cls._metrics_file) as metrics_file:
                db = yaml.load(metrics_file.read(), Loader=yaml.SafeLoader)
                db['tasks'].append(record)
            with open(cls._metrics_file, 'w') as metrics_file:
                metrics_file.write(yaml.dump(db))


class CloudCustomCoasterBackend(BaseBackend):
    def __init__(self, hostname='localhost', port=6000, runtask_py_path='runtask.py', **kwargs):
        super().__init__()
        self._task_id_counter = 0
        
        self.port = port
        self.hostname = hostname
        self.runtask_py_path = runtask_py_path
        
        self.total_cpus_in_pool = 800  # TODO make this dynamic
        self.bip = 10  # TODO make this synamic
        self._cpu_sem = BoundedSemaphore(int(self.total_cpus_in_pool))
        
        self.cloud_storage = AzureStorageSession(
            **kwargs['azure_params'],
            create_temp_container=True
            # create_temp_container=False
        )
        
        self.project_dir = os.getcwd()
        
        # Make directory for cjs launch scripts
        self.cloud_scripts_dir = os.path.join(self.project_dir, 'cloud_scripts')
        os.makedirs(self.cloud_scripts_dir, exist_ok=True)
        
        # Make directory for cjs launch scripts
        self.cloud_logs_dir = os.path.join(self.project_dir, 'cloud_logs')
        os.makedirs(self.cloud_logs_dir, exist_ok=True)
        
        # Initialize cloud metrics log
        CloudMetricsLogger.init(at=os.path.join(self.project_dir, 'cloud_metrics_{}.yaml'.format(datetime.now().strftime('%Y%m%d%H%M%S'))))
        
        
    
    def new_task_id(self, prefix='t'):
        self._task_id_counter += 1
        return prefix + str(self._task_id_counter).zfill(4)
    
    async def spawn(self, task):
        # Ensure the command body exists, otherwise there is nothing to do
        if 'cmd' not in task.directives:
            return task.complete()

        # Ensure there will ever be enough CPUs to run this task, otherwise fail
        task_requested_cpus = task.directives.get('cpus', 0)
        if task_requested_cpus > self.total_cpus_in_pool:
            log.critical('Task requested cpus ({}) greater than total available cpus ({})'.format(
                task_requested_cpus, self.total_cpus_in_pool
            ))
            return task.fail(1)
        
        # Determine whether this task should be run locally or on a remote cloud worker
        is_local_task = task.directives.get('cloud_args', dict()).get('local_task', False)
        log.info('Spawn ({}): {}'.format('Local' if is_local_task else 'Cloud', task))
        
        if is_local_task:
            return await self.spawn_local(task)
        
        # This is a cloud task
        return await self.spawn_cloud(task)
    
    async def spawn_cloud(self, task):
        task_id = self.new_task_id()
        start_time = datetime.now()
        bytes_sent_bundle, bytes_received_bundle = list(), list()
        
        try:
            fd_paths = {
                fd_name: fd
                for fd_name, fd in zip(('stdin', 'stdout', 'stderr'), self.get_fd_paths(task))
            }
            fd_filepointers = {
                fd_name: open(fd, fd_mode) if fd else None
                for fd_mode, (fd_name, fd) in zip(('r', 'w', 'w'), fd_paths.items())
            }
        except Exception as e:
            log.error(f'Exception in connecting input/output streams for task {task.name}: {e}')
            log.error(traceback.format_exc())
            fd_paths = {'stdin': None, 'stdout': None, 'stderr': None}
            fd_filepointers = {'stdin': None, 'stdout': None, 'stderr': None}
        
        try:
            data_metrics, ref_metrics, container_metrics = dict(), dict(), dict()
            data_metrics = blob_inputs_to_remote(task.directives['input'], self.cloud_storage)
            ref_inputs = parse_reference_input(task.directives.get('cloud_args', dict()).get('reference_input', list()))
            ref_metrics = blob_inputs_to_remote(ref_inputs, self.cloud_storage, blob_basename=True)  
            singularity_container_uri = task.directives.get('cloud_args', dict()).get('singularity_container')
            container_input = (
                [singularity_container_uri]
                if (
                    singularity_container_uri is not None
                    and not is_remote_uri(singularity_container_uri)
                    and get_cloud_directive('transfer_container_to_remote', task.directives)
                )
                else list()
            )
            container_metrics = blob_inputs_to_remote(container_input, self.cloud_storage)
            
            # Log metrics for input data
            total_input_metrics = data_metrics + ref_metrics + container_metrics
            for m in total_input_metrics:
                bytes_sent_bundle.append(m)
            bytes_sent_bundle.append({
                'name': 'total',
                'size': sum([max(0, t['size']) for t in total_input_metrics]),
                'time': sum([max(0, t['time']) for t in total_input_metrics])
            })
        except Exception as e:
            log.error(f'Exception in uploading data to the cloud bucket for task {task.name}: {e}')
            log.error(traceback.format_exc())
            
        
        try:
            remote_cmd = construct_customcoaster_cmd(
                task_body=task.directives['cmd'],
                cloud_downloads=task.directives['input'] + ref_inputs + container_input,
                cloud_uploads=task.directives['output'],
                account_name=self.cloud_storage.storage_account_name,
                account_key=self.cloud_storage.storage_account_key,
                cloud_scripts_dir=self.cloud_scripts_dir,
                container=self.cloud_storage._temp_container_name,
                singularity_container_uri=singularity_container_uri,
                task_name=task.name,
                runtask_py_path=self.runtask_py_path
            )
        except Exception as e:
            log.error(f'Cound not form remote worker command for task {task.name}: {e}')
            log.error(traceback.format_exc())
            remote_cmd = None
            
        try:
            p = await self.subprocess_sh(remote_cmd, **fd_filepointers)
            
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
                output_metrics = blob_outputs_to_local(task.directives['output'], self.cloud_storage)
                log.info(f'Complete: {task.name}')
                
                # Log metrics for output data
                for m in output_metrics:
                    bytes_received_bundle.append(m)
                bytes_received_bundle.append({
                    'name': 'total',
                    'size': sum([max(0, t['size']) for t in output_metrics]),
                    'time': sum([max(0, t['time']) for t in output_metrics])
                })
                
                return task.complete(p.returncode)
        except CancelledError:
            task.state['err'] = 'Runner cancelled Backend.spawn()'
            return task.fail(-15)
        finally:
            for fp in fd_filepointers.values():
                if fp is not None:
                    fp.close()
            
            # Get task runtime and which node it ran on
            elapsed_time = datetime.now() - start_time
            try:
                hostname = None
                with open(f'.{task.name}.hostname', 'r') as hostname_log:
                    hostname = hostname_log.read().strip()
                subprocess.call(['mv'] +  glob.glob('*.remote.out') + [self.cloud_logs_dir])
                subprocess.call(['mv'] +  glob.glob('*.remote.err') + [self.cloud_logs_dir])
                os.remove(f'.{task.name}.hostname')
            except:
                pass  # Fail silently
            
            CloudMetricsLogger.write_record({
                'task': task.name,
                'start_datetime': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'elapsed_time': str(elapsed_time),
                'end_datetime': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'in_files': bytes_sent_bundle,
                'out_files': bytes_received_bundle,
                'node': hostname or 'UNKNOWN'
            })

        return task
        
        # # Write out task body
        # with open(f'task{task_id}.sh', 'w') as task_body_sh:
        #     task_body_sh.write(task.directives['cmd'])
        # 
        # # Write out file list
        # # TODO Should this include the container?
        # with open(f'task{task_id}.files', 'w') as task_files:
        #     task_inputs = task.directives['input']
        #     ref_inputs = parse_reference_input(task.directives.get('cloud_args', dict()).get('reference_input', list()))
        #     for task_input in task_inputs + ref_inputs:
        #         task_files.write('\t'.join([
        #             task_input,
        #             'in',
        #             '',
        #             '',
        #             '',
        #             os.stat(task_input).st_size
        #         ]) + '\n')
        # 
        #     task_outputs = task.directives['output']
        #     for task_output in task_outputs:
        #         task_files.write('\t'.join([
        #             task_output,
        #             'out',
        #             '',
        #             '',
        #             '',
        #             os.stat(task_input).st_size
        #         ]) + '\n')
            
            
            
            
        
        
    async def spawn_local(self, task):
        cpus_reserved = 0
        fd_filepointers = dict()

        try:
            for i in range(task.directives.get('cpus', 0)):
                await self._cpu_sem.acquire()
                cpus_reserved += 1

            fd_paths = {
                fd_name: fd
                for fd_name, fd in zip(('stdin', 'stdout', 'stderr'), self.get_fd_paths(task))
            }
            fd_filepointers = {
                fd_name: open(fd, fd_mode) if fd else None
                for fd_mode, (fd_name, fd) in zip(('r', 'w', 'w'), fd_paths.items())
            }

            p = await self.subprocess_sh(
                task.directives['cmd'],
                **fd_filepointers
            )

            task.state.update(
                stdout_path=fd_paths['stdout'],
                stderr_path=fd_paths['stderr'],
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
            for fp in fd_filepointers.values():
                if fp is not None:
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