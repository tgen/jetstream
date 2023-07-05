import asyncio
import itertools
import json
import logging
import os
import sys
import re
import shlex
import shutil
import subprocess
import tempfile
import time
from asyncio import Lock
from asyncio.subprocess import PIPE
from datetime import datetime, timedelta
from jetstream.backends import BaseBackend
from jetstream.tasks import get_fd_paths
from jetstream import settings

log = logging.getLogger('jetstream.slurm')
SLURM_SACCT_DELIMITER = '\037'
SLURM_JOB_ID_PATTERN = re.compile(r"^(?P<jobid>\d+)(_(?P<arraystepid>\d+))?(\.(?P<stepid>(\d+|batch|extern)))?$")
SLURM_ACTIVE_STATES = settings['slurm_active_states'].get(list)
SLURM_PASSED_STATES = settings['slurm_passed_states'].get(list)
SLURM_SBATCH_RETRY = settings['slurm_sbatch_retry'].get(int)


class SlurmBackend(BaseBackend):
    """SlurmBackend will spawn tasks using a Slurm batch scheduler.

    The spawn coroutine will return when the slurm job ID for a task is
    complete. This works by maintaining a dict of SlurmBatchJobs, and
    periodically asking for updates from sacct."""
    count = itertools.count()
    respects = ('cmd', 'stdin', 'stdout', 'stderr', 'cpus', 'mem', 'walltime',
                'slurm_args')

    def __init__(
            self,
            sacct_frequency=60,
            sbatch_args=None,
            sbatch_delay=0.1,
            sbatch_executable=None,
            sacct_fields=('JobID', 'Elapsed', 'State', 'ExitCode'),
            job_monitor_max_fails=5):
        """SlurmBackend submits tasks as jobs to a Slurm batch cluster

        :param sacct_frequency: Frequency in seconds that job updates will
        be requested from sacct
        :param sbatch: path to the sbatch binary if not on PATH
        """
        super(SlurmBackend, self).__init__()
        self.sbatch_args = sbatch_args
        self.sbatch_executable = sbatch_executable
        self.sacct_frequency = sacct_frequency
        self.sacct_fields = sacct_fields
        self.sbatch_delay = sbatch_delay
        self.sbatch_lock = Lock()
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

        log.info('SlurmBackend initialized')

    def _bump_next_update(self):
        self._next_update = datetime.now() + timedelta(seconds=self.sacct_frequency)
        log.debug(f'Next sacct update bumped to {self._next_update.isoformat()}')

    def _get_sbatch_args(self, task):
        """Any extra args for sbatch will come from the application
        settings, followed by task settings. This means task settings
        will be able to override application config settings."""
        sbatch_args = []

        conf_sbatch_args = self.sbatch_args
        if conf_sbatch_args is None:
            pass
        elif isinstance(conf_sbatch_args, str):\
            # It can be a pain to store args s
            args = shlex.split(conf_sbatch_args)
            sbatch_args.extend(args)
        else:
            sbatch_args.extend(conf_sbatch_args)

        task_sbatch_args = task.directives.get('sbatch_args')
        if task_sbatch_args is None:
            pass
        elif isinstance(task_sbatch_args, str):
            args = shlex.split(task_sbatch_args)
            sbatch_args.extend(args)
        else:
            sbatch_args.extend(task_sbatch_args)

        return sbatch_args

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
                    sacct_data = sacct(*self.jobs, sacct_fields=self.sacct_fields, return_data=True)
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

        # sbatch breaks when called too frequently, so this places
        # a hard limit on the frequency of sbatch calls.
        time.sleep(self.sbatch_delay)

        stdin, stdout, stderr = get_fd_paths(task, self.runner.project)

        async with self.sbatch_lock:
            time.sleep(self.sbatch_delay)
            job = sbatch(
                cmd=task.directives['cmd'],
                name=task.name,
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
        log.info(f'SlurmBackend submitted({job.jid}): {task.name}')

        if sys.version_info > (3, 8):
            job.event = asyncio.Event()
        else:
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

            if state not in SLURM_ACTIVE_STATES:
                return True

        return False

    def is_ok(self):
        if not self.is_done():
            raise ValueError('Job is not complete yet.')

        if self.job_data['State'] in SLURM_PASSED_STATES:
            return True
        else:
            return False


def wait(*job_ids, sacct_fields=None, update_frequency=10):
    """Wait for one or more slurm batch jobs to complete"""
    while 1:
        jobs = sacct(*job_ids, sacct_fields=sacct_fields)

        if all([j.is_done() for j in jobs]):
            return
        else:
            time.sleep(update_frequency)


def sacct(*job_ids, sacct_fields=None, chunk_size=1000, strict=False, return_data=False):
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
        sacct_output = launch_sacct(*chunk, sacct_fields=sacct_fields)
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


def launch_sacct(*job_ids, sacct_fields=None, delimiter=SLURM_SACCT_DELIMITER, raw=False):
    """Launch sacct command and return stdout data

    This function returns raw query results, sacct() will be more
    useful in many cases.

    :param job_ids: Job ids to include in the query
    :param delimiter: Delimiter to separate parsable results data
    :param raw: Return raw stdout instead of parsed
    :return: Dict or Bytes
    """
    log.debug('Sacct request for {} jobs...'.format(len(job_ids)))
    args = ['sacct', '-P', '--format', '{}'.format(','.join(sacct_fields)), '--delimiter={}'.format(delimiter)]

    for jid in job_ids:
        args.extend(['-j', str(jid)])

    log.debug('Launching: {}'.format(' '.join([shlex.quote(r) for r in args])))
    p = subprocess.run(args, stdout=PIPE, check=True)

    if raw:
        return p.stdout.decode()

    return parse_sacct(p.stdout.decode(), delimiter=delimiter)


def parse_sacct(data, delimiter=SLURM_SACCT_DELIMITER, id_pattern=SLURM_JOB_ID_PATTERN):
    """Parse stdout from sacct to a dictionary of job ids and data."""
    jobs = dict()
    lines = iter(data.splitlines())
    header = next(lines).split(delimiter)

    for line in lines:
        row = dict(zip(header, line.split(delimiter)))

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


def sbatch(cmd, name=None, stdin=None, stdout=None, stderr=None, tasks=None,
           cpus_per_task=None, mem=None, walltime=None, comment=None,
           additional_args=None, sbatch_executable=None):
    if sbatch_executable is None:
        sbatch_executable = 'sbatch'

    args = [sbatch_executable, '--parsable']

    if name:
        args.extend(['-J', name])

    if stdin:
        args.extend(['--input', stdin])

    if stdout:
        args.extend(['-o', stdout])

    if stderr:
        args.extend(['-e', stdout])

    if tasks:
        args.extend(['-n', tasks])

    if cpus_per_task:
        args.extend(['-c', cpus_per_task])

    if mem:
        args.extend(['--mem', mem])

    if walltime:
        args.extend(['-t', walltime])

    if comment:
        args.extend(['--comment', comment])

    if additional_args:
        if isinstance(additional_args, str):
            args.append(additional_args)
        else:
            args.extend(additional_args)

    if cmd.startswith('#!'):
        script = cmd
    else:
        script = '#!/bin/bash\n{}'.format(cmd)

    temp = tempfile.NamedTemporaryFile()
    with open(temp.name, 'w') as fp:
        fp.write(script)

    args.append(temp.name)
    args = [str(r) for r in args]
    remaining_tries = SLURM_SBATCH_RETRY

    while 1:
        try:
            p = subprocess.run(args, stdout=subprocess.PIPE, check=True)
            break
        except subprocess.CalledProcessError:
            if remaining_tries == 0:
                raise
            else: 
                remaining_tries -= 1
                log.exception(f'Error during sbatch, retrying in 60s ...')
                time.sleep(60)

    jid = p.stdout.decode().strip()
    job = SlurmBatchJob(jid)
    job.args = args
    job.script = script
    return job
