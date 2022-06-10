import asyncio
import itertools
import json
import logging
import re
import shlex
import shutil
import subprocess
import tempfile
import time
from asyncio.subprocess import PIPE
from datetime import datetime, timedelta
from jetstream.backends import BaseBackend
from jetstream.tasks import get_fd_paths
from jetstream import settings

import dxpy
import dxpy.scripts
import dxpy.scripts.dx
from dxpy.utils.job_log_client import DXJobLogStreamClient

log = logging.getLogger('jetstream.backends')

def get_dx_instance( task ):
    name = task.directives.get("name")
    cpus_required = task.directives.get("cpus")
    memory_gb_required_str = task.directives.get('mem', "4G")
    memory_gb_required_value = int( memory_gb_required_str[:-1] )
    memory_gb_required_unit = memory_gb_required_str[-1]
    if memory_gb_required_unit == "M":
        memory_gb_required_value = 1
    elif memory_gb_required_unit != "G":
        raise RuntimeError('Task memory units must be M or G')
    memory_gb_required_per_cpu = float(memory_gb_required_value) / float( cpus_required )

    if memory_gb_required_per_cpu <= 2.0:
        mem_tier = "mem1"
    elif memory_gb_required_per_cpu <= 4.0:
        mem_tier = "mem2"
    else:
        mem_tier = "mem3"

    mem_cpu_tiers = { "mem1" : ( 2, 4, 8, 16, 36, 72 ),
                      "mem2" : ( 2, 4, 8, 16, 32, 48, 64, 96 ),
                      "mem3" : ( 2, 4, 8, 16, 32, 48, 64, 96 ) }

    for cpu_tier in mem_cpu_tiers[ mem_tier ]:
        if cpu_tier >= cpus_required:
            break
    
    cpu_tier = f'x{cpu_tier}'

    return f"{mem_tier}_ssd1_v2_{cpu_tier}"


class DnanexusBackend(BaseBackend):
    count = itertools.count()
    respects = ('cmd', 'stdin', 'stdout', 'stderr', 'cpus', 'mem', 'walltime',
                'dnanexus_args')

    def __init__(self, sacct_frequency=0.1, sbatch_delay=0.1,
                 sbatch_executable=None, sacct_fields=('JobID', 'Elapsed'),
                 job_monitor_max_fails=5,
                 dnanexus_app_name = 'docker-runner' ):
        super(DnanexusBackend, self).__init__()
        self.sbatch_executable = sbatch_executable
        self.sacct_frequency = sacct_frequency
        self.sacct_fields = sacct_fields
        self.sbatch_delay = sbatch_delay
        self.job_monitor_max_fails = job_monitor_max_fails
        self.jobs = dict()

        self.coroutines = (self.job_monitor,)
        self._next_update = datetime.now()

        if self.sbatch_executable is None:
            self.sbatch_executable = shutil.which('sbatch') or 'sbatch'

        self.dnanexus_app_name = "docker-runner"
        
        log.info('DnanexusBackend initialized')

    def _bump_next_update(self):
        self._next_update = datetime.now() + timedelta(seconds=self.sacct_frequency)
        log.debug(f'Next sacct update bumped to {self._next_update.isoformat()}')

    async def wait_for_next_update(self):
        """This allows the wait time to be bumped up each time a job is
        submitted. This means that sacct will never be checked immediately
        after submitting jobs, and it protects against finding data from
        old jobs with the same job ID in the database"""
        while datetime.now() < self._next_update:
            sleep_delta = self._next_update - datetime.now()
            sleep_seconds = max(0, sleep_delta.total_seconds())
            await asyncio.sleep(sleep_seconds)

    async def job_monitor(self):
        """Request job data updates from sacct for each job in self.jobs."""
        log.info('Dnanexus job monitor: Started!')
        failures = self.job_monitor_max_fails
        try:
            while 1:
                log.debug(f'Dnanexus job monitor: Waiting for next update')
                await self.wait_for_next_update()

                if not self.jobs:
                    log.debug('Dnanexus job monitor: No current jobs to check, bumping next update')
                    self._bump_next_update()
                    continue
                try:
                    log.debug(f'Dnanexus job monitor: Running sacct')
                    sacct_data = sacct(*self.jobs, return_data=True)
                except Exception:
                    if failures <= 0:
                        raise
                    else:
                        failures -= 1
                    err = f'Dnanexus job monitor: ERROR: {failures} remaining'
                    log.exception(err)
                    await asyncio.sleep(3)
                    continue

                failures = self.job_monitor_max_fails
                self._bump_next_update()

                for jid, data in sacct_data.items():
                    if jid in self.jobs:
                        job = self.jobs[jid]
                        job.job_data = data

                        if job.is_done() and job.event is not None:
                            log.debug(f'Dnanexus job monitor: Job done: {jid}')
                            job.event.set()
                            self.jobs.pop(jid)
        finally:
            log.info('Dnanexus job monitor: Stopped!')

    def dnanexus_job_name(self, task):
        """The dnanexus backend gives each job a name that is
        <run_id>.<job number>
        """
        count = next(self.count)
        run_id = self.runner.run_id
        return '{}.{}'.format(run_id, count)

    def dnanexus_job_comment(self, task):
        """Dnanexus jobs will receive a comment that contains details about the
        task, run id, and tags taken from the task directives. If tags are a
        string, they will be converted to a list with shlex.split"""
        tags = task.directives.get('tags', [])
        if isinstance(tags, str):
            tags = shlex.split(tags)

        comment = {
            'id': task.identity,
            'tags': tags
        }

        comment_string = json.dumps(comment, sort_keys=True)
        return comment_string

    def cancel(self):
        if self.jobs:
            log.info(f'Requesting scancel for {len(self.jobs)} dnanexus jobs')
            for job_id, job in self.jobs.items():
                job.dx_job.terminate()
           
    async def spawn(self, task):
        log.debug(f'Spawn: {task.name}')

        if not task.directives.get('cmd'):
            return task.complete()

        time.sleep(self.sbatch_delay)
        stdin, stdout, stderr = get_fd_paths(task, self.runner.project)

        container = task.directives.get( 'container' )
        instance_type = get_dx_instance( task )
        
        job = sbatch(
            cmd=task.directives['cmd'],
            container=container,
            instance_type=instance_type,
            input_filenames=task.directives.get( 'input', [] ),
            output_filenames=task.directives.get( 'output', [] ),
        )

        task.state.update(
            label = f'Dnanexus({job.jid})',
            stdout_path = stdout,
            stderr_path = stderr,
            dnanexus_job_id = job.jid,
            dnanexus_cmd = job.args
        )

        self._bump_next_update()
        log.info(f'DnanexusBackend submitted({job.jid}): {task.name}')

        job.event = asyncio.Event()
        self.jobs[ job.jid ] = job

        await job.event.wait()
        log.debug('Job event was set, gathering info updating status')

        if self.sacct_fields:
            job_info = {k: v for k, v in job.job_data.items() if
                        k in self.sacct_fields}
            task.state['dnanexus_sacct'] = job_info
        
        if stdout is not None:
            with open(stdout, 'w') as outfile:
              log_client = DXJobLogStreamClient( job.jid,
                                                 print_job_info = False,
                                                 msg_callback = lambda m: outfile.write( m.get( "msg" ) )  )
              log_client.connect()
    
        if job.is_ok():
            log.info(f'Complete: {task.name}')
            task.complete(job.returncode())
        else:
            log.info(f'Failed: {task.name}')
            task.fail(job.returncode())

        log.debug(f'DnanexusBackend spawn completed for {task.name}')
        return task


class DnanexusBatchJob(object):
    states = {
        'idle': 'Job is idle',
        'runnable': 'Job is runnable',
        'running': 'Job is running',
        'done': 'Job is done',
        'failed': 'Job has failed'
    }

    active_states = {'idle', 'runnable', 'running'}

    inactive_states = {'done', 'failed'}

    failed_states = {'failed'}

    passed_states = {'done'}

    def __init__(self, dx_job_id):
        self.dx_job = dxpy.bindings.dxjob.DXJob( dxid = dx_job_id )
        dx_job_describe = self.dx_job.describe()
        self.args = [ dx_job_describe["runInput"]["command_string"] ]
        self.script = dx_job_describe["runInput"]["command_string"]
        self.jid = dx_job_describe["id"]
        
        self._job_data = None

    def __eq__(self, other):
        try:
            if other == self.jid:
                return True
        except AttributeError:
            return other.jid == self.jid

    def __repr__(self):
        return '<DnanexusBatchJob: {}>'.format(self.jid)

    def _update_state(self, job_data):
        self._job_data = job_data

    def update(self):
        data = launch_sacct(self.jid)

        if not self.jid in data:
            raise ValueError('No job data found for:  {}'.format(self.jid))
        else:
            self.job_data = data[self.jid]

    def wait(self, *args, **kwargs):
        return wait(self.jid, *args, **kwargs)

    @property
    def job_data(self):
        return self._job_data

    @job_data.setter
    def job_data(self, value):
        self._update_state(value)

    def returncode(self):
        if not self.is_done():
            raise ValueError('Job not done yet')
        try:
            return int(self.job_data['ExitCode'].partition(':')[0])
        except (KeyError, IndexError):
            if self.is_ok():
                return 0
            else:
                return 1

    def cancel(self):
        log.info( f'Terminating {self.jid}')
        self.dx_job.terminate()
        return "OK"

    def is_done(self):
        if self._job_data:
            state = self._job_data.get('State')

            if state not in self.active_states:
                return True

        return False

    def is_ok(self):
        if not self.is_done():
            raise ValueError('Job is not complete yet.')

        if self._job_data['State'] in self.passed_states:
            return True
        else:
            return False


def wait(*job_ids, update_frequency=10):
    """Wait for one or more dnanexus batch jobs to complete"""
    while 1:
        jobs = sacct(*job_ids)

        if all([j.is_done() for j in jobs]):
            return
        else:
            time.sleep(update_frequency)


def sacct(*job_ids, chunk_size=1000, strict=False, return_data=False):
    if not job_ids:
        raise ValueError('Missing required argument: job_ids')

    job_ids = [str(jid) for jid in job_ids]
    jobs = [DnanexusBatchJob(jid) for jid in job_ids]

    data = {}
    for i in range(0, len(job_ids), chunk_size):
        chunk = job_ids[i: i + chunk_size]
        sacct_output = launch_sacct(*chunk)
        data.update(sacct_output)

    log.debug('Status updates for {} jobs'.format(len(data)))

    if return_data:
        return data

    for job in jobs:
        if not job.jid in data:
            if strict:
                raise ValueError('No records returned for {}'.format(job.jid))
            else:
                log.debug('No records found for {}'.format(job.jid))
        else:
            job.job_data = data[job.jid]

    return jobs


def launch_sacct(*job_ids, raw=False):
    """Launch sacct command and return stdout data

    This function returns raw query results, sacct() will be more
    useful in many cases.

    :param job_ids: Job ids to include in the query
    :param delimiter: Delimiter to separate parsable results data
    :param raw: Return raw stdout instead of parsed
    :return: Dict or Bytes
    """
    log.debug('Sacct request for {} jobs...'.format(len(job_ids)))
    
    jobs = {}
    for jid in job_ids:
        dx_job = dxpy.bindings.dxjob.DXJob( dxid = jid )
        dx_job_describe = dx_job.describe()
        dx_job_id = dx_job_describe["id"]
        state = dx_job_describe["state"]
        
        log.debug( f'JOB {jid} {dx_job_id} {state}' )
        jobs[ dx_job_id ] = { "State" : state }
    
    return jobs

def normalize_command_string( cmd ):
    cmd_lines = []
    for i, line in enumerate( cmd.splitlines() ):
        line = re.sub( r'\\\\', 'JETSTREAM_DOUBLE_BACKSLASH', line )
        line = re.sub( '#.*', '', line )
        cmd_lines.append( line )
        
    normalized_command_string = "; ".join( cmd_lines )
    normalized_command_string = re.sub( r'\\;(\s*)', '\g<1>', normalized_command_string )
    normalized_command_string = re.sub( 'JETSTREAM_DOUBLE_BACKSLASH', r'\\\\', normalized_command_string )
    normalized_command_string = re.sub( '(;\s*)+', '; ', normalized_command_string )
    normalized_command_string = re.sub( '^\s*;', '', normalized_command_string )
    return normalized_command_string

def sbatch( cmd,
            container,
            instance_type,
            input_filenames = [],
            output_filenames = [],
            dnanexus_project_id = None,
            dnanexus_project_work_directory = None
          ):
  
    dnanexus_project_id = None
    dnanexus_project_work_directory = None
    dnanexus_app_name = 'docker-runner'
    
    if dnanexus_project_id is None:
        workspace_describe = dxpy.describe( dxpy.WORKSPACE_ID )
        dnanexus_project_id = workspace_describe[ 'id' ]
    
    if dnanexus_project_work_directory is None:
        dnanexus_project_work_directory = dxpy.scripts.dx.get_pwd().split(":")[1]
    
    # command_string = normalize_command_string( cmd )
    command_string = cmd
    
    dx_app_run_args = { 'command_string' : command_string,
                        'docker_image' : container,
                        'dnanexus_project_id' : dnanexus_project_id,
                        'dnanexus_project_work_directory' : dnanexus_project_work_directory,
                        'input_filenames' : input_filenames,
                        'output_filenames' : output_filenames }

    dx_app = dxpy.bindings.dxapp.DXApp( name = dnanexus_app_name )

    try:
        dx_job = dx_app.run( dx_app_run_args, instance_type = instance_type  )
    except Exception as e:
        log.warning(f'''\
Failed to run dx job:
command string: {dx_app_run_args.command_string}
docker image: {dx_app_run_args.docker_image}
''')
        return None

    dx_job_id = dx_job.describe()["id"]
    
    job = DnanexusBatchJob( dx_job_id )
    return job


if __name__ == "__main__":
    
    dnanexus_project_id = None
    dnanexus_project_work_directory = None
    dnanexus_app_name = 'docker-runner'
    
    if dnanexus_project_id is None:
        workspace_describe = dxpy.describe( dxpy.WORKSPACE_ID )
        dnanexus_project_id = workspace_describe[ 'id' ]
    
    if dnanexus_project_work_directory is None:
        dnanexus_project_work_directory = dxpy.scripts.dx.get_pwd().split(":")[1]
        
    
    new_dx_file_args = { "project" : dnanexus_project_id,
                         "folder" : dnanexus_project_work_directory,
                         "name" : "test.json" }
    
    my_new_dx_file = dxpy.DXFile()
    my_new_dx_file.new(**new_dx_file_args)
    dxpy.upload_local_file( filename = "test.json",
                            wait_on_close = True,
                            use_existing_dxfile = my_new_dx_file )
    my_new_dx_file.close( block=True )
