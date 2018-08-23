import os
import re
import shlex
import time
import json
import subprocess
from datetime import datetime, timedelta
import itertools
from multiprocessing import cpu_count
import asyncio
from asyncio import BoundedSemaphore, Event, create_subprocess_shell
from asyncio.subprocess import PIPE
from concurrent.futures import CancelledError
from jetstream.workflows import save_workflow
from jetstream import utils
from jetstream import log


class AsyncRunner(object):
    def __init__(self, backend=None, max_forks=None, output_prefix='',
                 autosave=0, default_yield=1):
        """AsyncRunner executes a workflow using a given backend
        :param backend: Backend object used to spawn tasks
        :param output_prefix:
        :param max_forks: If this is omitted, an estimate of 25% of the system
            thread limit is used. Over subscribing the fork limit may cause the
            system to reject creation of new processes, and tasks will fail to
            launch.
        """
        self.backend = backend or LocalBackend()
        self.backend.runner = self
        self.output_prefix = output_prefix
        self.autosave = timedelta(seconds=autosave)
        self.default_yield = default_yield
        self.workflow = None
        self.project = None
        self.fp = None
        self.run = Event()
        self._start_time = None
        self._loop = None
        self._workflow_manager = None
        self._last_save = datetime.now()
        self._sem = BoundedSemaphore(max_forks or guess_max_forks())

    def log_status(self):
        """ Logs a status report """
        log.info('AsyncRunner event-loop load: {}'.format(
            len(asyncio.Task.all_tasks())))
        log.info('Workflow status: {}'.format(self.workflow))

    async def _yield(self, delay=None):
        if delay is None:
            delay = self.default_yield

        log.verbose('Yield for {}s'.format(delay))

        if self.autosave and (datetime.now() - self._last_save) > self.autosave:
            self._last_save = datetime.now()

            if self.project:
                path = self.project.workflow_file
            else:
                # TODO move this to a config module
                path = 'jetstream_workflow_{}.yaml'.format(self.fp.id)

            save_workflow(self.workflow, path=path)

        await asyncio.sleep(delay)

    async def workflow_manager(self):
        """Workflow manager will constantly ask workflow for new tasks
        and spawn a new coroutine for each task when ready. """
        log.info('Workflow manager started!')

        try:
            for task in self.workflow:
                if task is None:
                    await self._yield()
                else:
                    await self.spawn(task)

                await self._yield(0)
        except Exception as e:
            log.exception(e)
            log.info('Exception in workflow manager, halting run!')
            self._loop.stop()
        finally:
            log.info('Workflow manager stopped!')
            self.run.clear()

    async def spawn(self, task):
        log.verbose('Registering backend.spawn: {}'.format(task))
        log.verbose('Runner concurrency semaphore: {}'.format(self._sem))

        if task.get('cmd') is None:
            return task.complete(0)

        try:
            await self._sem.acquire()
            asyncio.ensure_future(self.backend.spawn(task))
        except Exception as e:
            log.exception(e)
            log.info('Exception during task spawn, halting run!')
            self._loop.stop()
        finally:
            self._sem.release()

    def start_backend_coros(self):
        if hasattr(self.backend, 'start_coros'):
            coros = self._loop.create_task(self.backend.start_coros())
        else:
            coros = self._loop.create_task(asyncio.sleep(0))

        return coros

    def finalize_run(self):
        fails = [t for t in self.workflow.tasks(objs=True) if
                 t.status == 'failed']

        if fails:
            log.info('\u2620  Some tasks failed! {}'.format(fails))
            rc = 1
        else:
            rc = 0

        return rc

    def close(self):
        if self._loop:
            self._loop.close()

    def start(self, workflow, project=None, loop=None):
        self.workflow = workflow
        self.project = project
        self.fp = utils.Fingerprint()
        self.run.set()

        start = datetime.now()
        os.environ.update(JETSTREAM_RUN_ID=self.fp.id)
        log.info('Starting Run ID: {}'.format(self.fp.id))

        if project:
            history_file = os.path.join(project.history_dir, self.fp.id)
            with open(history_file, 'w') as fp:
                utils.yaml_dump(self.fp.serialize(), stream=fp)

        if loop is None:
            self._loop = asyncio.get_event_loop()
        else:
            self._loop = loop

        try:
            manager = self._loop.create_task(self.workflow_manager())
            coros = self.start_backend_coros()
            self._loop.run_until_complete(manager)
            coros.cancel()
        except KeyboardInterrupt:
            log.critical('Received interrupt: Shutting down!')
            self.run.clear()
            for t in asyncio.Task.all_tasks():
                t.cancel()
        finally:
            rc = self.finalize_run()
            log.info('Shutting down event loop')
            self._loop.run_until_complete(self._loop.shutdown_asyncgens())

        elapsed = datetime.now() - start
        log.info('Run {} Elapsed: {}'.format(self.fp.id, elapsed))
        return rc


class Backend(object):
    """Backends should implement the coroutine method "spawn".

    Backend.spawn will be called by the AsyncRunner when a task is ready
    for execution. It should be asynchronous and return the task_id
    and return code as a tuple: (task_id, returncode)

    Backends should declare a set of task directives that are used when
    executing a task: Backend.respects """
    respects = set()
    runner = None

    async def create_subprocess_shell(self, cmd, **kwargs):
        while self.runner.run.is_set():
            try:
                return await create_subprocess_shell(cmd, **kwargs)
            except BlockingIOError as e:
                log.info('Unable to start subprocess: {}'.format(e))
                log.info('Retry in 10 seconds.')
                await asyncio.sleep(10)

    async def subprocess_run_sh(
            self, args, *, stdin=None, input=None, stdout=None, stderr=None,
            cwd=None, check=False, encoding=None, errors=None, env=None,
            loop=None, executable='/bin/bash'):
        """Asynchronous version of subprocess.run

        This will always use a shell to launch the subprocess, and it prefers
        /bin/bash (can be changed via arguments)"""
        log.verbose('subprocess_run_sh: {}'.format(args))

        p = await self.create_subprocess_shell(
            args, stdin=stdin, stdout=stdout, stderr=stderr, cwd=cwd,
            encoding=encoding, errors=errors, env=env,
            loop=loop, executable=executable)

        stdout, stderr = await p.communicate(input=input)

        if check and p.returncode != 0:
            raise ChildProcessError(args)

        return subprocess.CompletedProcess(args, p.returncode, stdout, stderr)


class LocalBackend(Backend):
    respects = ('cmd', 'stdin', 'stdout', 'stderr', 'cpus')
    
    def __init__(self, cpus=None):
        """The LocalBackend executes tasks as processes on the local machine.
        
        :param cpus: If this is None, the number of available CPUs will be
            guessed.
        """
        n_cpus = cpus or guess_local_cpus()
        self._cpu_sem = BoundedSemaphore(n_cpus)
        log.info('LocalBackend initialized with {} cpus'.format(n_cpus))
    
    def status(self):
        return 'CPUs: {}'.format(self._cpu_sem)

    async def spawn(self, task):
        log.info('Spawn: {}'.format(task))
        log.debug(self.status())

        cpus_reserved = 0

        try:
            for i in range(task.get('cpus', 0)):
                await self._cpu_sem.acquire()
                cpus_reserved += 1

            if task.get('stdin'):
                input = task['stdin'].encode()
            else:
                input = None

            cmd = task['cmd']
            
            if task.get('stdout'):
                stdout = self.runner.output_prefix + task.get('stdout')
            else:
                stdout = None
                
            if task.get('stderr'):
                stderr = self.runner.output_prefix + task.get('stderr')
            else:
                stderr = None
            
            if stdout:
                stdout_fp = open(stdout, 'w')
            else:
                stdout_fp = None
            
            if stderr:
                stderr_fp = open(stderr, 'w')
            else:
                stderr_fp = None
            
            p = await self.subprocess_run_sh(
                cmd, stdin=PIPE, input=input, stdout=stdout_fp, stderr=stderr_fp
            )

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
        log.info('Sacct request for {} jobs...'.format(len(self._jobs)))
        sacct_data = {}

        for chunk in self.chunk_jobs():
            data = await self.async_sacct_request(*chunk)
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

    async def async_sacct_request(self, *job_ids):
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

    def sacct_request(self, *job_ids):
        if not job_ids:
            raise ValueError('Missing required argument "job_ids"')

        job_args = ' '.join(['-j {}'.format(jid) for jid in job_ids])

        cmd = 'sacct -P --format all --delimiter={} {}'.format(
            self.sacct_delimiter, job_args)

        log.debug('Launching: {}'.format(cmd))
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
        run_id = self.runner.fp['id']
        count = next(self.count)
        job_name = '{}.{}'.format(run_id, count)

        comment = json.dumps({
            'run': self.runner.fp.serialize(),
            'task': task.serialize()
        }, sort_keys=True)

        args = ['sbatch', '--parsable', '-J', job_name, '--comment', comment]
        
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

            log.info("Submitted batch job {}".format(jid))

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
                # self.returncode = self.job_data['ExitCode'].partition(':')[0]

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
