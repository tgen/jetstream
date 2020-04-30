import asyncio
import glob
import itertools
import json
import logging
import os
import pathlib
import re
import shlex
import shutil
import subprocess
import tempfile
import time
from asyncio import Lock, BoundedSemaphore, create_subprocess_shell, CancelledError
from asyncio.subprocess import PIPE
from datetime import datetime, timedelta
from jetstream.backends import BaseBackend
from jetstream import settings

log = logging.getLogger('jetstream.slurm')
sacct_delimiter = '\037'
job_id_pattern = re.compile(r"^(?P<jobid>\d+)(_(?P<arraystepid>\d+))?(\.(?P<stepid>(\d+|batch|extern)))?$")

class SlurmSingularityBackend(BaseBackend):
    """SlurmSingularityBackend will spawn tasks using a Slurm batch scheduler.

    The spawn coroutine will return when the slurm job ID for a task is
    complete. This works by maintaining a dict of SlurmBatchJobs, and
    periodically asking for updates from sacct."""
    count = itertools.count()
    respects = ('cmd', 'stdin', 'stdout', 'stderr', 'cpus', 'mem', 'walltime',
                'slurm_args')

    def __init__(self, sacct_frequency=5, sbatch_delay=0.1,
                 sbatch_executable=None, sacct_fields=('JobID', 'Elapsed'),
                 job_monitor_max_fails=5):
        """SlurmSingularityBackend submits tasks as jobs to a Slurm batch cluster

        :param sacct_frequency: Frequency in seconds that job updates will
        be requested from sacct
        :param sbatch: path to the sbatch binary if not on PATH
        """
        super(SlurmSingularityBackend, self).__init__()
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

        with open(os.devnull, 'w') as devnull:
            subprocess.run(
                [self.sbatch_executable, '--version'],
                check=True,
                stdout=devnull,
                stderr=devnull
            )

        # self._singularity_run_sem = BoundedSemaphore( self.cpus )
        self._singularity_pull_lock = Lock()
        self._singularity_pull_cache = {}
        
        log.info('SlurmSingularityBackend initialized')

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
        log.info('Slurm job monitor started!')
        failures = self.job_monitor_max_fails
        try:
            while 1:
                await self.wait_for_next_update()

                if not self.jobs:
                    log.debug('No current jobs to check')
                    self._bump_next_update()
                    continue
                try:
                    sacct_data = sacct(*self.jobs, return_data=True)
                except Exception:
                    if failures <= 0:
                        raise
                    else:
                        failures -= 1
                    err = f'Slurm job monitor error: {failures} remaining'
                    log.exception(err)
                    await asyncio.sleep(120)
                    continue

                failures = self.job_monitor_max_fails
                self._bump_next_update()

                for jid, data in sacct_data.items():
                    if jid in self.jobs:
                        job = self.jobs[jid]
                        job.job_data = data

                        if job.is_done() and job.event is not None:
                            job.event.set()
                            self.jobs.pop(jid)
        finally:
            log.info('Slurm job monitor stopped!')

    def slurm_job_name(self, task):
        """The slurm backend gives each job a name that is
        <run_id>.<job number>
        """
        count = next(self.count)
        run_id = self.runner.run_id
        return '{}.{}'.format(run_id, count)

    def slurm_job_comment(self, task):
        """Slurm jobs will receive a comment that contains details about the
        task, run id, and tags taken from the task directives. If tags are a
        string, they will be converted to a list with shlex.split"""
        tags = task.directives.get('tags', [])
        if isinstance(tags, str):
            tags = shlex.split(tags)

        comment = {
            'id': task.identity,
            'tags': tags
        }

        # TODO long tags here could cause an sbatch submission error, but
        # we still have over 1000 characters before that happens, so the
        # chance is pretty small. Limiting the tag length is non-trivial
        # but possible. It would probably be best to enforce this limit at
        # the Task level though.

        comment_string = json.dumps(comment, sort_keys=True)
        return comment_string

    def cancel(self):
        if self.jobs:
            jobs = list(self.jobs.keys())
            log.info(f'Requesting scancel for {len(jobs)} slurm jobs')
            subprocess.run(['scancel'] + jobs)

    async def spawn(self, task):
        log.debug(f'Spawn: {task.name}')

        if not task.directives.get('cmd'):
            return task.complete()

        # sbatch breaks when called too frequently
        time.sleep(self.sbatch_delay)
        stdin, stdout, stderr = self.get_fd_paths(task)

        input_filenames = task.directives.get( 'input', [] )
        output_filenames = task.directives.get( 'output', [] )
        
        docker_image = task.directives.get( 'docker_image', None )
        if docker_image == None:
            raise RuntimeError(f'docker_image argument missing for task: {task.name}')
        
        singularity_image = f"docker://{docker_image}"
        
        log.debug( f'going to pull: {singularity_image}' )
        try:
            if singularity_image in self._singularity_pull_cache:
                pass
            else:
                async with self._singularity_pull_lock:
                    # for i in range( self.cpus ):
                    #     await self._singularity_run_sem.acquire()
                    opt_https = "--nohttps " if singularity_image.startswith("docker://localhost") else ""
                    pull_command_run_string = f'singularity exec {opt_https}{singularity_image} true'
                    log.debug( f'pulling: {pull_command_run_string}' )
                    _p = await create_subprocess_shell( pull_command_run_string,
                                                        stdout=asyncio.subprocess.PIPE,
                                                        stderr=asyncio.subprocess.PIPE )
                    stdout, stderr = await _p.communicate()
                    log.debug( f'pulled, stdout: {stdout}' )
                    log.debug( f'pulled, stderr: {stderr}' )
                    self._singularity_pull_cache[ singularity_image ] = singularity_image
                    # for i in range( self.cpus ):
                    #     self._singularity_run_sem.release()
        except Exception as e:
            log.warning(f'Exception during singularity prepull: {e}')
            p = await create_subprocess_shell( "exit 1;" )


        log.debug( f'pulled: {singularity_image}' )
        
        job = sbatch(
            cmd=task.directives['cmd'],
            singularity_image=singularity_image,
            name=task.name,
            input_filenames=input_filenames,
            output_filenames=output_filenames,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            comment=self.slurm_job_comment(task),
            cpus_per_task=task.directives.get('cpus'),
            mem=task.directives.get('mem'),
            walltime=task.directives.get('walltime'),
            additional_args=task.directives.get('sbatch_args'),
            sbatch_executable=self.sbatch_executable
        )

        task.state.update(
            label=f'Slurm({job.jid})',
            stdout_path=stdout,
            stderr_path=stderr,
            slurm_job_id=job.jid,
            slurm_cmd=' '.join(shlex.quote(a) for a in job.args)
        )

        self._bump_next_update()
        log.info(f'SlurmSingularityBackend submitted({job.jid}): {task.name}')

        job.event = asyncio.Event(loop=self.runner.loop)
        self.jobs[job.jid] = job

        await job.event.wait()
        log.debug(f'{task.name}: job info was updated')

        if self.sacct_fields:
            job_info = {k: v for k, v in job.job_data.items() if
                        k in self.sacct_fields}
            task.state['slurm_sacct'] = job_info

        if job.is_ok():
            log.info(f'Complete: {task.name}')
            task.complete(job.returncode())
        else:
            log.info(f'Failed: {task.name}')
            task.fail(job.returncode())

        log.debug(f'Slurmbackend returning task: {task.name}')
        return task


class SlurmBatchJob(object):
    states = {
        'BOOT_FAIL': 'Job terminated due to launch failure, typically due to a '
                     'hardware failure (e.g. unable to boot the node or block '
                     'and the job can not be requeued).',
        'CANCELLED': 'Job was explicitly cancelled by the user or system '
                     'administrator. The job may or may not have been '
                     'initiated.',
        'COMPLETED': 'Job has terminated all processes on all nodes with an '
                     'exit code of zero.',
        'CONFIGURING': 'Job has been allocated resources, but are waiting for '
                       'them to become ready for use (e.g. booting).',
        'COMPLETING': 'Job is in the process of completing. Some processes on '
                      'some nodes may still be active.',
        'FAILED': 'Job terminated with non-zero exit code or other failure '
                  'condition.',
        'NODE_FAIL': 'Job terminated due to failure of one or more allocated '
                     'nodes.',
        'PENDING': 'Job is awaiting resource allocation.',
        'PREEMPTED': 'Job terminated due to preemption.',
        'REVOKED': 'Sibling was removed from cluster due to other cluster '
                   'starting the job.',
        'RUNNING': 'Job currently has an allocation.',
        'SPECIAL_EXIT': 'The job was requeued in a special state. This state '
                        'can be set by users, typically in EpilogSlurmctld, if '
                        'the job has terminated with a particular exit value.',
        'STOPPED': 'Job has an allocation, but execution has been stopped with '
                   'SIGSTOP signal. CPUS have been retained by this job.',
        'SUSPENDED': 'Job has an allocation, but execution has been suspended '
                     'and CPUs have been released for other jobs.',
        'TIMEOUT': 'Job terminated upon reaching its time limit.'
    }

    active_states = {'CONFIGURING', 'COMPLETING', 'RUNNING', 'SPECIAL_EXIT',
                     'PENDING'}

    inactive_states = {'BOOT_FAIL', 'CANCELLED', 'COMPLETED', 'FAILED',
                       'NODE_FAIL', 'PREEMPTED', 'REVOKED',
                       'STOPPED', 'SUSPENDED', 'TIMEOUT'}

    failed_states = {'BOOT_FAIL', 'CANCELLED', 'FAILED', 'NODE_FAIL'}

    passed_states = {'COMPLETED'}

    def __init__(self, jid=None, data=None):
        self.args = None
        self._job_data = None

        if data:
            if jid is None:
                self.jid = str(data['JobID'])
            self._update_state(data)
        else:
            self.jid = str(jid)

    def __eq__(self, other):
        try:
            if other == self.jid:
                return True
        except AttributeError:
            return other.jid == self.jid

    def __repr__(self):
        return '<SlurmBatchJob: {}>'.format(self.jid)

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
        """Attempts to returns a standard integer exit code based on Slurm
        "derived" exit code, but falls back to some dumb heuristics if the
        exit code isn't parsing."""
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
        log.info('Launching "scancel {}"'.format(self.jid))
        cmd_args = ('scancel', self.jid)
        return subprocess.call(cmd_args)

    def is_done(self):
        if self._job_data:
            state = self._job_data.get('State')

            if state not in self.active_states:
                return True

        return False

    def is_ok(self):
        if not self.is_done():
            raise ValueError('Job is not complete yet.')

        if self.job_data['State'] in self.passed_states:
            return True
        else:
            return False


def wait(*job_ids, update_frequency=10):
    """Wait for one or more slurm batch jobs to complete"""
    while 1:
        jobs = sacct(*job_ids)

        if all([j.is_done() for j in jobs]):
            return
        else:
            time.sleep(update_frequency)


def sacct(*job_ids, chunk_size=1000, strict=False, return_data=False):
    """Query sacct for job records.

    Jobs are returned for each job id, but steps will be combined under a
    single job id object. This will return a placeholder job for any job
    id given, regardless of whether job data was returned by sacct. The
    strict option can be used to raise an error when job data is missing
    for any of the job ids."""
    if not job_ids:
        raise ValueError('Missing required argument: job_ids')

    job_ids = [str(jid) for jid in job_ids]
    jobs = [SlurmBatchJob(jid) for jid in job_ids]

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


def launch_sacct(*job_ids, delimiter=sacct_delimiter, raw=False):
    """Launch sacct command and return stdout data

    This function returns raw query results, sacct() will be more
    useful in many cases.

    :param job_ids: Job ids to include in the query
    :param delimiter: Delimiter to separate parsable results data
    :param raw: Return raw stdout instead of parsed
    :return: Dict or Bytes
    """
    log.debug('Sacct request for {} jobs...'.format(len(job_ids)))
    args = ['sacct', '-P', '--format', 'all', '--delimiter={}'.format(delimiter)]

    for jid in job_ids:
        args.extend(['-j', str(jid)])

    log.debug('Launching: {}'.format(' '.join([shlex.quote(r) for r in args])))
    p = subprocess.run(args, stdout=PIPE, check=True)

    if raw:
        return p.stdout.decode()

    return parse_sacct(p.stdout.decode(), delimiter=delimiter)


def parse_sacct(data, delimiter=sacct_delimiter, id_pattern=job_id_pattern):
    """Parse stdout from sacct to a dictionary of job ids and data."""
    jobs = dict()
    lines = iter(data.strip().splitlines())
    header = next(lines).strip().split(delimiter)

    for line in lines:
        row = dict(zip(header, line.strip().split(delimiter)))

        try:
            match = id_pattern.match(row['JobID'])
            groups = match.groupdict()
        except (KeyError, AttributeError):
            # Job data restrictions are very loose - there is a small chance
            # that the chosen delimiter was added to some field in the job
            # data, and that will break this parser. If that happens records
            # are skipped and a warning is issued, but parsing continues.
            log.warning('Error parsing sacct line: {}'.format(line))
            continue

        # Slurm job ids are <jobid>[_<arrayid>][.<taskid>]. The goal here
        # is to group all job steps (tasks, array steps) under their
        # corresponding jid. The steps are added to a list under the key
        # "_steps", all other data updates the dictionary.
        jid = groups['jobid']

        if groups['stepid'] or groups['arraystepid']:
            if jid not in jobs:
                jobs[jid] = {'_steps': list()}

            jobs[jid]['_steps'].append(row)
        else:
            if jid not in jobs:
                jobs[jid] = dict()

            row['_steps'] = list()
            jobs[jid].update(row)

    return jobs


def sbatch(cmd, singularity_image, name=None, input_filenames=[], output_filenames=[],
           stdin=None, stdout=None, stderr=None, tasks=None,
           cpus_per_task=1, mem="2 G", walltime="1h", comment=None,
           additional_args=None, sbatch_executable=None, retry=10):
    
    # determine input/output mounts needed
    singularity_mounts = set()
    for input_filename_glob_pattern in input_filenames:
        input_filenames_glob = glob.glob( input_filename_glob_pattern )
        if len( input_filenames_glob ) == 0:
            raise RuntimeError(f'input file(s) do not exist: {input_filename_glob_pattern}')
        for input_filename in input_filenames_glob:
            input_filename = os.path.abspath( input_filename )
            input_filename_head, input_filename_tail = os.path.split( input_filename )
            singularity_mounts.add( input_filename_head )
    for output_filename in output_filenames:
        output_filename = os.path.abspath( output_filename )
        output_filename_head, output_filename_tail = os.path.split( output_filename )
        singularity_mounts.add( output_filename_head )
        pathlib.Path( output_filename_head ).mkdir( parents=True, exist_ok=True )
    mount_strings = []
    for singularity_mount in singularity_mounts:
        mount_strings.append( "-B %s" % ( singularity_mount ) )
    singularity_mounts_string = " ".join( mount_strings )
    
    # create cmd script
    os.makedirs( "jetstream/cmd", mode = 0o777, exist_ok = True )
    millis = int(round(time.time() * 1000))
    if name == None:
        name = "script"
    cmd_script_filename = f"jetstream/cmd/{millis}_{name}.cmd"
    cmd_script_filename = os.path.abspath( cmd_script_filename )
    with open( cmd_script_filename, "w" ) as cmd_script:
        cmd_script.write( cmd )
    
    # create sbatch script
    # if sbatch_executable is None:
    #     sbatch_executable = 'sbatch'
    # 
    # sbatch_args = [sbatch_executable, '--parsable']
    sbatch_args = []

    if name:
        sbatch_args.extend(['-J', name])

    if stdin:
        sbatch_args.extend(['--input', stdin])

    if stdout:
        sbatch_args.extend(['-o', stdout])

    if stderr:
        sbatch_args.extend(['-e', stdout])

    if tasks:
        sbatch_args.extend(['-n', tasks])

    if cpus_per_task:
        sbatch_args.extend(['-c', cpus_per_task])

    if mem:
        sbatch_args.extend(['--mem', mem])

    if walltime:
        sbatch_args.extend(['-t', walltime])

    if comment:
        sbatch_args.extend(['--comment', comment])

    if additional_args:
        if isinstance(additional_args, str):
            sbatch_args.append(additional_args)
        else:
            sbatch_args.extend(additional_args)
    
    sbatch_script = "#!/bin/bash\n"
    for i in range( 0, len(sbatch_args), 2 ):
        sbatch_script += f"#SBATCH {sbatch_args[i]} {sbatch_args[i+1]}\n"
        
    opt_https = "--nohttps " if singularity_image.startswith("docker://localhost") else ""
    sbatch_script += f"#!/bin/bash\nsingularity exec --nv {opt_https}{singularity_mounts_string} {singularity_image} bash {cmd_script_filename}\n"

    if name == None:
        name = "script"
    sbatch_script_filename = f"jetstream/cmd/{millis}_{name}.sbatch"
    sbatch_script_filename = os.path.abspath( sbatch_script_filename )
    with open( sbatch_script_filename, "w" ) as sbatch_script_file:
        sbatch_script_file.write( sbatch_script )
        
    # command_run_string = f"""singularity exec --nv {opt_https}{singularity_mounts_string} {singularity_image} bash {run_script_filename}"""
    # script = '#!/bin/bash\n{}'.format(command_run_string)

    # create sbatch script
    # opt_https = "--nohttps " if singularity_image.startswith("docker://localhost") else ""
    # run_singularity_script = """#!/bin/bash\nsingularity exec --nv {opt_https}{singularity_mounts_string} {singularity_image} bash {cmd_script_filename}"""
    # 
    # temp = tempfile.NamedTemporaryFile()
    # with open(temp.name, 'w') as fp:
    #     fp.write(run_singularity_script)
    #     
    # args.append(command_run_string)
    # args = [str(r) for r in args]
    
    submit_sbatch_args = [ "sbatch", f"{sbatch_script_filename}" ]
    
    remaining_tries = int(retry)
    while 1:
        try:
            p = subprocess.run(submit_sbatch_args, stdout=subprocess.PIPE, check=True)
            break
        except subprocess.CalledProcessError:
            if remaining_tries > 0:
                remaining_tries -= 1
                log.exception(f'Error during sbatch, retrying in 60s ...')
                time.sleep(60)
            else:
                raise

    jid = p.stdout.decode().strip().split()[-1]
    job = SlurmBatchJob(jid)
    job.args = submit_sbatch_args
    job.script = sbatch_script_filename
    return job
