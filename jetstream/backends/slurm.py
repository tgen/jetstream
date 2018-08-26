import re
import shlex
import time
import json
import subprocess
import itertools
import asyncio
from asyncio import BoundedSemaphore, Event
from asyncio.subprocess import PIPE
from concurrent.futures import CancelledError
from jetstream import log
from jetstream.backends import Backend


class SlurmBackend(Backend):
    count = itertools.count()
    sacct_delimiter = '\037'
    submission_pattern = re.compile(r"Submitted batch job (\d*)")
    job_id_pattern = re.compile(r"(?P<jobid>\d*)\.?(?P<taskid>.*)")
    respects = ('cmd', 'stdin', 'stdout', 'stderr', 'cpus', 'mem', 'walltime',
                'slurm_args')

    def __init__(self, max_jobs=9001, sacct_frequency=10, chunk_size=1000):
        self.sacct_frequency = sacct_frequency
        self.chunk_size = chunk_size
        self._jobs = {}
        self._jobs_sem = BoundedSemaphore(max_jobs)

        log.info('SlurmBackend initialized with {} max jobs'.format(max_jobs))

    def status(self):
        return 'Slurm jobs: {}'.format(self._jobs_sem)

    def chunk_jobs(self):
        seq = list(self._jobs.keys())
        size = self.chunk_size
        return (seq[pos:pos + size] for pos in range(0, len(seq), size))

    def add_jid(self, jid):
        job = SlurmBatchJob(jid)
        self._jobs[jid] = job
        return job

    async def start_coros(self):
        log.info('Slurm job monitor started!')

        try:
            while self.runner.run.is_set():
                await asyncio.sleep(self.sacct_frequency)
                await self._update_jobs()
        except CancelledError:
            if self._jobs:
                jids = ' '.join(list(self._jobs.keys()))
                log.info('Requesting scancel for: {}'.format(jids))
                await self.subprocess_run_sh('scancel {}'.format(jids))
        finally:
            log.info('Slurm job monitor stopped!')

    async def _update_jobs(self):
        log.verbose('Sacct request for {} jobs...'.format(len(self._jobs)))
        sacct_data = {}

        for chunk in self.chunk_jobs():
            data = await self.async_sacct_request(*chunk)
            sacct_data.update(data)

        log.verbose('Status updates for {} jobs'.format(len(sacct_data)))

        reap = set()
        for jid, job in self._jobs.items():
            if jid in sacct_data:
                log.debug('Updating: {}'.format(jid))
                job_data = sacct_data[jid]
                job.update(job_data)

                if job.is_complete:
                    reap.add(jid)
            else:
                log.warning('No sacct data found for {}'.format(jid))

        for jid in reap:
            self._jobs.pop(jid)

    async def async_sacct_request(self, *job_ids):
        if not job_ids:
            raise ValueError('Missing required argument "job_ids"')

        job_args = ' '.join(['-j {}'.format(jid) for jid in job_ids])

        cmd = 'sacct -P --format all --delimiter={} {}'.format(
            self.sacct_delimiter, job_args)

        log.verbose('Launching: {}'.format(cmd))
        p = await self.create_subprocess_shell(cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = await p.communicate()

        res = self._parse_sacct(stdout.decode())

        return res

    def sacct_request(self, *job_ids):
        if not job_ids:
            raise ValueError('Missing required argument "job_ids"')

        job_args = ' '.join(['-j {}'.format(jid) for jid in job_ids])

        cmd = 'sacct -P --format all --delimiter={} {}'.format(
            self.sacct_delimiter, job_args)

        log.verbose('Launching: {}'.format(cmd))
        stdout = subprocess.check_output(cmd)

        res = self._parse_sacct(stdout.decode())

        return res

    def _parse_sacct(self, data):
        """ Parse stdout from sacct to a dictionary of jobs and job_data. """
        if not data:
            return {}

        jobs = dict()
        lines = iter(data.splitlines())
        header = next(lines).strip().split(self.sacct_delimiter)

        for line in lines:
            row = dict(zip(header, line.strip().split(self.sacct_delimiter)))
            match = self.job_id_pattern.match(row['JobID'])

            if match is None:
                log.info('Unable to parse sacct line: {}'.format(line))
                pass

            groups = match.groupdict()
            jobid = groups['jobid']
            taskid = groups['taskid']

            if taskid is '':
                if jobid in jobs:
                    log.info('Duplicate record for job: {}'.format(jobid))
                else:
                    row['_steps'] = list()
                    jobs[jobid] = row
            else:
                if jobid not in jobs:
                    jobs[jobid] = {'_steps': list()}

                jobs[jobid]['_steps'].append(row)

        log.debug('Parsed data for {} jobs'.format(len(jobs)))
        return jobs

    def sbatch_cmd(self, task):
        """ Returns a formatted sbatch command. """
        run_id = self.runner.fp.id
        count = next(self.count)
        job_name = '{}.{}'.format(run_id, count)

        tags = task.get('tags', [])
        if isinstance(tags, str):
            tags = tags.split()

        comment = json.dumps({
            'run': self.runner.fp.serialize(),
            'task': {
                'tid': task.tid,
                'tags': tags,
            }
        }, sort_keys=True)

        args = ['sbatch', '--parsable', '-J', job_name,
                '--comment \'{}\''.format(comment)]

        if task.get('stdout'):
            if task.get('stderr'):
                stdout = self.runner.output_prefix + task['stdout']
                stderr = self.runner.output_prefix + task['stderr']
                args.extend(['-o', stdout, '-e', stderr])
            else:
                stdout = self.runner.output_prefix + task['stdout']
                args.extend(['-o', stdout])
        else:
            if task.get('stderr'):
                stderr = self.runner.output_prefix + task['stderr']
                args.extend(['-e', stderr])
            else:
                pass  # Don't add any o/e args to slurm command

        # Slurm requires that we request at least 1 cpu
        if task.get('cpus'):
            args.extend(['-c', str(task['cpus'])])

        if task.get('mem'):
            args.extend(['--mem', str(task['mem'])])

        if task.get('walltime'):
            args.extend(['-t', str(task['walltime'])])

        if task.get('slurm_args'):
            args.extend(task['slurm_args'])

        if task.get('stdin'):
            formatted = 'echo \'{}\' | {}'.format(task['stdin'], task['cmd'])
            final_cmd = shlex.quote(formatted)
        else:
            final_cmd = shlex.quote(task['cmd'])

        args.extend(['--wrap', final_cmd])

        return args

    async def spawn(self, task):
        log.debug('Spawn: {}'.format(task))
        job = None

        try:
            await self._jobs_sem.acquire()

            # Sbatch fails when called too frequently so here is a bandaid
            time.sleep(.1)

            sbatch_cmd = self.sbatch_cmd(task)
            cmd = ' '.join(sbatch_cmd)

            log.debug('Final command: {}'.format(cmd))

            p = await self.subprocess_run_sh(cmd, stdout=PIPE)

            jid = p.stdout.decode().strip()

            if p.returncode != 0:
                log.info('Error submitting job: {}'.format(jid))
                return task.fail(1)

            log.info("{} Slurm JobID: {}".format(task, jid))
            task.set_state(slurm_job_id=jid)

            job = self.add_jid(jid)
            rc = await job.wait()

            if rc != 0:
                return task.fail(rc)
            else:
                return task.complete(rc)

        except CancelledError:
            if job is not None:
                await self.subprocess_run_sh('scancel {}'.format(job.jid))

            return task.complete(-15)

        finally:
            self._jobs_sem.release()


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

    def __init__(self, jid):
        self.jid = str(jid)
        self._job_data = None
        self._returncode = None
        self._is_complete = Event()

    @property
    def is_complete(self):
        return self._is_complete.is_set()

    @property
    def job_data(self):
        return self._job_data

    @job_data.setter
    def job_data(self, value):
        self._job_data = value
        state = self._job_data.get('State', '')

        if state not in self.active_states:
            try:
                if state not in self.passed_states:
                    self.returncode = 1
                else:
                    self.returncode = 0

                    # TODO
                    # Slurm sets the return code to 0 when it cancels jobs for
                    # memory issues etc.. I consider this a failure, but this
                    # may need to be revisited later.

            except KeyError:
                self.returncode = -123

    def update(self, job_data):
        self.job_data = job_data

    @property
    def returncode(self):
        return self._returncode

    @returncode.setter
    def returncode(self, value):
        self._returncode = int(value)
        self._is_complete.set()

    async def wait(self):
        try:
            await self._is_complete.wait()
        except CancelledError:
            self.cancel()

        return self.returncode

    def cancel(self):
        log.info('Scancel: {}'.format(self.jid))
        cmd_args = ('scancel', self.jid)
        return subprocess.call(cmd_args)
