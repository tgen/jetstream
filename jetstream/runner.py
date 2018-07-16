import os
import re
import logging
import shlex
import time
import subprocess
import jetstream
from datetime import datetime
from multiprocessing import cpu_count
import asyncio
from asyncio import (BoundedSemaphore, Event, create_subprocess_shell,
                     create_subprocess_exec)
from asyncio.subprocess import PIPE, STDOUT
from concurrent.futures import CancelledError

log = logging.getLogger(__name__)


class AsyncRunner(object):
    def __init__(self, workflow, backend=None, logging_interval=None,
                 log_path=None):
        self.workflow = workflow
        self.backend = backend
        self.backend.runner = self
        self.logging_interval = logging_interval or 3
        self.log_path = log_path or 'logs'
        self._fp = jetstream.utils.Fingerprint()

        self._started = False
        self._loop = None
        self._iterator = None
        self._workflow_manager = None
        self._logger = None
        self._workflow_complete = Event()

    @property
    def fp(self):
        return self._fp

    def log_status(self):
        """ Logs a status report """
        log.info('AsyncRunner event-loop load: {}'.format(
            len(asyncio.Task.all_tasks())))
        log.info('Workflow status: {}'.format(self.workflow))

    async def logger(self):
        log.info('Logger started!')

        try:
            while not self._workflow_complete.is_set():
                await asyncio.sleep(self.logging_interval)
                self.log_status()
        finally:
            log.info('Logger stopped!')

    async def workflow_manager(self):
        log.info('Workflow manager started!')

        self._iterator = iter(self.workflow)

        try:
            while 1:
                try:
                    task = next(self._iterator)
                except StopIteration:
                    break

                if task is None:
                    await asyncio.sleep(.1)
                else:
                    await self.spawn(task)

                await asyncio.sleep(0)
        finally:
            log.info('Workflow manager stopped!')
            self._workflow_complete.set()

    async def spawn(self, task):
        log.debug('Registering backend.spawn: {}'.format(task))

        try:
            asyncio.ensure_future(self.backend.spawn(task))
        except Exception as e:
            log.exception(e)
            log.info('Exception during task spawn, halting run!')
            self._loop.stop()

    def start_backend_coros(self):
        if hasattr(self.backend, 'start_coros'):
            coros = self._loop.create_task(self.backend.start_coros())
        else:
            coros = self._loop.create_task(asyncio.sleep(0))

        return coros

    def init_run(self):
        self._start_time = datetime.now()
        log.info('Run ID: {}'.format(self.fp.id))

    def finalize_run(self):
        self._loop.run_until_complete(self._loop.shutdown_asyncgens())
        self._loop.close()

        log.info('Finalizing run: {}'.format(self.fp.id))

        elapsed = datetime.now() - self._start_time
        log.info('Total run time: {}'.format(elapsed))


        fails = [t.id for t in self.workflow.tasks(objs=True) if
                 t.status == 'failed']

        if fails:
            log.info('\u2620  Some tasks failed! {}'.format(fails))
            rc = 1
        else:
            rc = 0

        self._start_time = None
        return rc

    def start(self, loop=None):
        if self._started:
            raise RuntimeError('Runners can only be started once')
        else:
            self._started = True

        try:
            self.init_run()

            if loop is None:
                self._loop = asyncio.get_event_loop()
            else:
                self._loop = loop

            manager = self._loop.create_task(self.workflow_manager())
            logger = self._loop.create_task(self.logger())
            coros = self.start_backend_coros()

            self._loop.run_until_complete(manager)
            logger.cancel()
            coros.cancel()
        except KeyboardInterrupt:
            for t in asyncio.Task.all_tasks():
                t.cancel()
        finally:
            return self.finalize_run()


class Backend(object):
    """ Backends should implement the coroutine method "spawn" """

    # "spawn" will be called by the asyncrunner when a task is ready
    # for execution. It should be asynchronous and return the task_id
    # and returncode as a tuple: (task_id, returncode)
    def __init__(self, max_forks=None):
        self.runner = None
        self._sem = BoundedSemaphore(max_forks or guess_max_forks())
        log.info('Initialized with: {}'.format(self._sem))

    def status(self):
        return 'Backend status: {}'.format(str(self._sem))

    @staticmethod
    async def create_subprocess_shell(cmd, **kwargs):
        while 1:
            try:
                return await create_subprocess_shell(cmd, **kwargs)
            except BlockingIOError as e:
                log.info('Unable to start subprocess: {}'.format(e))
                log.info('Retry in 10 seconds.')
                await asyncio.sleep(10)

    @staticmethod
    async def create_subprocess_exec(*args, **kwargs):
        while 1:
            try:
                return await create_subprocess_exec(*args, **kwargs)
            except BlockingIOError as e:
                log.info('Unable to start subprocess: {}'.format(e))
                log.info('Retry in 10 seconds.')
                await asyncio.sleep(10)

    async def subprocess_run(
            self, args, *, stdin=None, input=None, stdout=None, stderr=None,
            shell=False, cwd=None, check=False, encoding=None,
            errors=None, env=None, loop=None):
        """ Asynchronous version of subprocess.run """
        log.debug('Subprocess run: {}'.format(self._sem))

        await self._sem.acquire()

        try:
            if shell:
                p = await self.create_subprocess_shell(
                    args, stdin=stdin, stdout=stdout, stderr=stderr, cwd=cwd,
                    encoding=encoding, errors=errors, env=env, 
                    loop=loop)
            else:
                p = await self.create_subprocess_exec(
                    *args, stdin=stdin, stdout=stdout, stderr=stderr, cwd=cwd,
                    encoding=encoding, errors=errors, env=env, 
                    loop=loop)

            stdout, stderr = await p.communicate(input=input)

            if check and p.returncode != 0:
                raise ChildProcessError(args)

            c = subprocess.CompletedProcess(args, p.returncode, stdout, stderr)
            return c

        finally:
            self._sem.release()

    async def spawn(self, *args, **kwargs):
        raise NotImplementedError


class LocalBackend(Backend):
    def __init__(self, cpus=None, *args, **kwargs):
        """ LocalBackend executes tasks as processes on the local machine.
        
        :param max_subprocess: Total number of subprocess allowed regardless of
        cpu requirements. If this is omitted, an estimate of 25% of the system 
        thread limit is used. If child processes in a workflow spawn several 
        threads, the system may start to reject creation of new processes.
        """
        super(LocalBackend, self).__init__(*args, **kwargs)
        self._cpu_sem = BoundedSemaphore(cpus or guess_local_cpus())
        log.info('LocalBackend initialized with: {}'.format(self._cpu_sem))
    
    def status(self):
        return 'CPUs: {} Forks: {} '.format(self._cpu_sem, self._sem)

    async def spawn(self, task):
        log.info('Spawn: {}'.format(task))
        log.debug(self.status())

        cpus_reserved = 0

        if task.cmd is None:
            return task.complete(0)

        try:
            for i in range(task.cpus):
                await self._cpu_sem.acquire()
                cpus_reserved += 1

            if task.stdin:
                input = task.stdin.encode()
            else:
                input = None

            cmd = task.cmd
            stdout = os.path.join(self.runner.log_path, task.stdout)
            stderr = os.path.join(self.runner.log_path, task.stderr)

            if stdout == stderr:
                with open(stdout, 'w') as fp:
                    p = await self.subprocess_run(cmd, stdin=PIPE, input=input,
                                                  stdout=fp, stderr=STDOUT,
                                                  shell=True)
            else:
                with open(stdout, 'w') as out, open(stderr, 'w') as err:
                    p = await self.subprocess_run(cmd, stdin=PIPE, input=input,
                                                  stdout=out, stderr=err,
                                                  shell=True)

            if p.returncode != 0:
                return task.fail(p.returncode)
            else:
                return task.complete(p.returncode)

        except CancelledError:
            return task.fail(-15)

        finally:
            log.info('Done: {}'.format(task))

            for i in range(cpus_reserved):
                self._cpu_sem.release()


class SlurmBackend(Backend):
    sacct_delimiter = '\037'
    submission_pattern = re.compile(r"Submitted batch job (\d*)")
    job_id_pattern = re.compile(r"(?P<jobid>\d*)\.?(?P<taskid>.*)")

    def __init__(self, max_jobs=None, sacct_frequency=None, chunk_size=None,
                 *args, **kwargs):
        super(SlurmBackend, self).__init__(*args, **kwargs)
        self.sacct_frequency = sacct_frequency or 60
        self.chunk_size = chunk_size or 10000
        self._monitor_jobs = Event()
        self._jobs = {}
        self._jobs_sem = BoundedSemaphore(max_jobs or 1000000)

        log.info('SlurmBackend initialized with: {}'.format(self._jobs_sem))

    def status(self):
        return 'Forks: {} Slurm jobs: {}'.format(self._sem, self._jobs_sem)

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
        self._monitor_jobs.set()

        try:
            while self._monitor_jobs.is_set():
                await asyncio.sleep(self.sacct_frequency)
                await self._update_jobs()
        except CancelledError:
            self._monitor_jobs.clear()
        finally:
            log.info('Slurm job monitor stopped!')

    async def _update_jobs(self):
        log.info('Sacct request for {} jobs...'.format(len(self._jobs)))
        sacct_data = {}

        for chunk in self.chunk_jobs():
            data = await self._sacct_request(*chunk)
            sacct_data.update(data)

        log.info('Status updates for {} jobs'.format(len(sacct_data)))

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

    async def _sacct_request(self, *job_ids):
        if not job_ids:
            raise ValueError('Missing required argument "job_ids"')

        job_args = ' '.join(['-j {}'.format(jid) for jid in job_ids])

        cmd = 'sacct -P --format all --delimiter={} {}'.format(
            self.sacct_delimiter, job_args)

        log.debug('Launching: {}'.format(cmd))
        p = await self.create_subprocess_shell(cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = await p.communicate()

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
        args = ['sbatch', '--parsable', '-J', task.id, '--comment',
                shlex.quote(self.runner.fp.to_json())]

        stdout = os.path.join(self.runner.log_path, task.stdout)
        stderr = os.path.join(self.runner.log_path, task.stderr)

        if stdout == stderr:
            args.extend(['-o', stdout])
        else:
            args.extend(['-o', stdout, '-e', stderr])

        # Slurm requires that we request at least 1 cpu
        if task.cpus != 0:
            args.extend(['-c', str(task.cpus)])

        if task.mem != 0:
            args.extend(['--mem', str(task.mem)])

        if task.walltime != 0:
            args.extend(['-t', str(task.time)])

        if task.stdin:
            formatted_cmd = 'echo \'{}\' | {}'.format(task.stdin, task.cmd)
            final_cmd = shlex.quote(formatted_cmd)
        else:
            final_cmd = shlex.quote(task.cmd)

        args.extend(['--wrap', final_cmd])

        return args

    async def spawn(self, task):
        log.debug('Spawn: {}'.format(task))

        job = None

        try:
            await self._jobs_sem.acquire()

            # Sbatch fails when called too frequently so here is a bandaid
            time.sleep(.1)

            if task.cmd is None:
                return task.complete(0)

            sbatch_cmd = self.sbatch_cmd(task)
            cmd = ' '.join(sbatch_cmd)

            log.debug('Final command: {}'.format(cmd))

            p = await self.subprocess_run(cmd, stdout=PIPE, stderr=STDOUT,
                                          shell=True)

            jid = p.stdout.decode().strip()

            if p.returncode != 0:
                log.info('Error submitting job: {}'.format(jid))
                return task.fail(1)

            log.info("Submitted batch job {}".format(jid))

            job = self.add_jid(jid)
            rc = await job.wait()

            if rc != 0:
                return task.fail(rc)
            else:
                return task.complete(rc)

        except CancelledError:
            if job is not None:
                await self.subprocess_run('scancel {}'.format(job.jid))

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
                self.returncode = self.job_data['ExitCode'].partition(':')[0]
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


def guess_max_forks(default=500):
    try:
        res = int(0.25 * int(subprocess.check_output('ulimit -u', shell=True)))
        return res
    except FileNotFoundError as e:
        log.exception(e)
        return default


def guess_local_cpus(default=4):
    return cpu_count() or default
