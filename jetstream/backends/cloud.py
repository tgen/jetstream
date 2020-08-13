import logging
import asyncio
import jetstream
from jetstream.backends import BaseBackend
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
import yaml
import time
import traceback

from ..utils import construct_cjs_cmd, dynamic_import
from ..cloud.azure import AzureStorageSession
from ..cloud.base import blob_inputs_to_remote, blob_outputs_to_local, is_remote_uri, CLOUD_STORAGE_PROVIDERS

log = logging.getLogger('jetstream.cloud')

def get_pool_info(pool_name, api_key, pworks_url="http://beta.parallel.works"):
    """
    With a valid ParallelWorks API Key, get relevant information on the worker pool
    
    :param pool_name: str Name of an active pool in PW
    :param api_key: str Valid API for a PW account
    :param pworks_url: str URL of the PW platform
    :return: dict The requested pool information
    """
    rq = requests.get('{pworks_url}/api/resources?key={api_key}'.format(pworks_url=pworks_url, api_key=api_key))
    for pool_data in rq.json():
        if pool_data['name'] == pool_name:
            log.info('Got pooldata: {}'.format(pool_data))
            return {
                'serviceport': pool_data['info']['ports']['serviceport'],
                'controlport': pool_data['info']['ports']['controlport'],
                'maxworkers': int(pool_data['settings']['max']),
                # TODO Sometimes the pool information doesn't contain the info needed for this calculation, but 
                # it isn't really necessary to have anyway
                # 'cpus': int(pool_data['info']['cpuPerWorker']) // pool_data['settings']['jobsPerNode']
                'cpus': 8
            }


def parse_reference_input(ref_input_directive):
    """
    Allows for wildcards to be interpreted as blobs on the local filesystem
    
    :param ref_input_directive: list<str>|str Paths given in the cloud_args.ref_input directive
    :return: list<str> A list of paths with all wildcards expanded
    """
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
    """
    Helper function to get a directive one layer undernearth ``cloud_args``
    
    :param key: str Directive key
    :param task_directives: dict The dictionary of task directives
    :param cloud_args_key: str Key for the first level
    :return: object The requested ``cloud_args`` directive
    """
    return task_directives.get(cloud_args_key, dict()).get(key)


class CloudSwiftBackend(BaseBackend):
    """
    Executes tasks on cloud-based worker nodes, handling all data transfer to/from worker nodes and 
    the client node.
    """
    def __init__(self, pw_pool_name=None, pw_api_key=None, cpus=None, blocking_io_penalty=None,
                 cloud_storage_provider='azure', **kwargs):
        """
        If ``pw_pool_name`` and ``pw_api_key`` are both given valid values, this backend will use ParallelWorks to 
        manage resouces elastically. If they remain ``None``, then it is assumed that the use has manually started 
        a pool using the included ``start_pool.py`` script. Note that both approaches requires the use of binaries 
        available as part of the Swift workflow language.
        
        :param pw_pool_name: str Name of an active pool in PW
        :param pw_api_key: str Valid API for a PW account
        :param cpus: int The total number of CPUs available to the worker pool
        :param blocking_io_penalty: int Delay (in seconds) when a BlockingIOError prevents a new process from spawning.
        :param cloud_storage_provder: str Name of the cloud storage provider, which must match up with one of the keys 
            in ``jetstream.cloud.base.CLOUD_STORAGE_PROVIDERS``
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
        self.total_cpus_in_pool = self.pool_info['cpus'] * self.pool_info['maxworkers']
        self.bip = blocking_io_penalty \
                   or jetstream.settings['backends']['local']['blocking_io_penalty'].get(int)
        self._cpu_sem = BoundedSemaphore(int(self.total_cpus_in_pool))
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
        
        # Instantiate a cloud storage provider
        storage_class = dynamic_import(CLOUD_STORAGE_PROVIDERS[cloud_storage_provider])
        self.cloud_storage = storage_class(
            **kwargs[storage_class.config_key]
        )
        
        # Initialize cloud metrics log
        CloudMetricsLogger.init(at=os.path.join(self.project_dir, 'cloud_metrics_{}.yaml'.format(datetime.now().strftime('%Y%m%d%H%M%S'))))
        
        log.info(f'CloudSwiftBackend initialized with {self.total_cpus_in_pool} cpus')

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
        cmd = task.directives['cmd']
        cpus_reserved = 0
        
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
                cloud_storage=self.cloud_storage,
                cjs_stagein=None,
                cjs_stageout=None,
                cloud_downloads=task.directives['input'] + reference_inputs + container_input,
                cloud_uploads=task.directives['output'],
                cloud_scripts_dir=self.cloud_scripts_dir,
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


class CloudMetricsLogger:
    """
    Helper class to keep records on data transferred and tasks running on worker nodes.
    """
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
