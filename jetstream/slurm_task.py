""" Slurm tasks are jobs launched on a Slurm scheduler

Slurm-specific task directives:

mem
----

Specify the real memory required per node. Default units are megabytes unless
the SchedulerParameters configuration parameter includes the "default_gbytes"
option for gigabytes. Different units can be specified using the suffix
[K|M|G|T]. Default value is DefMemPerNode and the maximum value is
MaxMemPerNode. If configured, both parameters can be seen using the scontrol
show config command. This parameter would generally be used if whole nodes are
allocated to jobs (SelectType=select/linear). Also see --mem-per-cpu. --mem
and --mem-per-cpu are mutually exclusive.

NOTE: A memory size specification of zero is treated as a special case and
grants the job access to all of the memory on each node. If the job is
allocated multiple nodes in a heterogeneous cluster, the memory limit on each
node will be that of the node in the allocation with the smallest memory size
(same limit will apply to every node in the job's allocation).

NOTE: Enforcement of memory limits currently relies upon the task/cgroup
plugin or enabling of accounting, which samples memory use on a periodic basis
(data need not be stored, just collected). In both cases memory use is based
upon the job's Resident Set Size (RSS). A task may exceed the memory limit
until the next periodic accounting sample.

cpus
-----
Advise the Slurm controller that ensuing job steps will require ncpus number
of processors per task. Without this option, the controller will just try to
allocate one processor per task.

time
-----

Set a limit on the total run time of the job allocation. If the requested time
limit exceeds the partition's time limit, the job will be left in a PENDING
state (possibly indefinitely). The default time limit is the partition's
default time limit. When the time limit is reached, each task in each job step
is sent SIGTERM followed by SIGKILL. The interval between signals is specified
by the Slurm configuration parameter KillWait. The OverTimeLimit configuration
parameter may permit the job to run longer than scheduled. Time resolution is
one minute and second values are rounded up to the next minute.

A time limit of zero requests that no time limit be imposed. Acceptable time
formats include "minutes", "minutes:seconds", "hours:minutes:seconds",
"days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".


"""
import re
import time
import subprocess
import logging
from jetstream.runners import BaseTask

log = logging.getLogger(__name__)

slurm_max_req_freq = 30
slurm_sacct_delimiter = '\037'
slurm_submission_pattern = re.compile(r"Submitted batch job (\d*)")
slurm_job_id_pattern = re.compile(r"(?P<jobid>\d*)\.?(?P<taskid>.*)")
slurm_states = {
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


class SlurmTask(BaseTask):
    def __init__(self, task_id, task_directives):
        super(SlurmTask, self).__init__(task_id, task_directives)
        self._launched = False
        self._job = None

    def poll(self):
        return self._job.poll()

    def kill(self):
        return self._job.cancel()

    def wait(self):
        self.returncode = self._job.wait()
        return self.returncode

    def launch(self):
        self._launched = True
        args = ['-J', self.task_id]

        if self.stdout_path != self.stderr_path:
            args.extend(['-o', self.stdout_path, '-e', self.stderr_path])
        else:
            args.extend(['-o', self.stdout_path])

        if 'cpus' in self.task_directives:
            args.extend(['-c', self.task_directives['cpus']])

        if 'mem' in self.task_directives:
            args.extend(['--mem', self.task_directives['mem']])

        if 'time' in self.task_directives:
            args.extend(['-t', self.task_directives['time']])

        try:
            self._job = sbatch(*args, stdin_data=self.stdin_data)
            return True
        except Exception as e:
            log.exception(e)
            return False


class Sopen(object):
    _watching = set()
    _last_poll = 0

    def __init__(self, jid, cluster=None, job_data=None):
        """Tracks a Slurm job and provides methods similar to a Popen object"""
        Sopen._watching.add(self)

        self.jid = str(jid)
        self.cluster = cluster

        if job_data is not None:
            self.job_data = job_data
        else:
            self.job_data = dict()

    def __repr__(self):
        return "Sopen({})".format(self.jid)

    @property
    def status(self):
        try:
            state = self.job_data['State']
            if state not in slurm_states:
                return 'COMPLETED'
            else:
                return state
        except KeyError:
            return None

    @property
    def returncode(self):
        try:
            return int(self.job_data['ExitCode'].partition(':')[0])
        except KeyError:
            return None

    def poll(self):
        if self.status is None:
            self.update()

        if self.status in active_states:
            return None
        else:
            return self.returncode

    def cancel(self):
        return scancel('-j', self.jid)

    def wait(self, timeout=None):
        started = time.time()

        while 1:
            if self.status in inactive_states:
                break
            else:
                self.update()

                elapsed = time.time() - started
                if timeout and elapsed > timeout:
                    raise TimeoutError

        return self.returncode

    def update(self, force=False):
        """ Request a status update from Slurm. """
        if force:
            self.job_data = get_one(self.jid)

        elif Sopen._watching:
            # Requesting an update for any Sopen object will automatically
            # batch the update request along with any other Sopen instances.
            # This event is limited to one request every slurm_max_req_freq
            # seconds.

            now = time.time()
            if now - Sopen._last_poll < slurm_max_req_freq:
                return

            data = get_all(*[j.jid for j in Sopen._watching])

            _prune = []
            for job in Sopen._watching:
                log.critical('Updating: {}'.format(job))
                try:
                    job.job_data = data[job.jid]
                except KeyError:
                    _prune.append(job)

            for i in _prune:
                Sopen._watching.remove(i)

            if self.status in inactive_states:
                Sopen._watching.remove(self)


def _sacct_request(*job_ids):
    global slurm_last_poll

    delay_until = slurm_last_poll + slurm_max_req_freq
    time.sleep(max(0, delay_until - time.time()))

    if not job_ids:
        raise ValueError('Missing required argument "job_ids"')

    cmd = ['sacct', '-P', '--format', 'all',
           '--delimiter={}'.format(slurm_sacct_delimiter)]

    for jid in job_ids:
        cmd.extend(['-j', jid])

    log.critical('Launching: {}'.format(' '.join(cmd)))
    res = subprocess.check_output(cmd).decode()

    slurm_last_poll = time.time()
    return res


def _parse_sacct(data):
    jobs = dict()
    lines = iter(data.splitlines())
    header = next(lines).strip().split(slurm_sacct_delimiter)

    for line in lines:
        row = dict(zip(header, line.strip().split(slurm_sacct_delimiter)))
        m = slurm_job_id_pattern.match(row['JobID'])

        if m is None:
            log.critical('Unable to parse sacct line: {}'.format(line))
            pass

        groups = m.groupdict()
        jobid = groups['jobid']
        taskid = groups['taskid']

        if taskid is '':
            if jobid in jobs:
                log.critical('Duplicate record for job: {}'.format(jobid))
            else:
                row['_steps'] = list()
                jobs[jobid] = row
        else:
            if not jobid in jobs:
                jobs[jobid] = {'_steps': list()}

            jobs[jobid]['_steps'].append(row)

    log.critical('Parsed data for {} jobs.'.format(len(jobs)))
    return jobs


def get_one(job_id):
    """Get job data for a single job. This raises an exception if the
    job is not found in the sacct results. """
    job_id = str(job_id)
    data = _sacct_request(job_id)
    jobs = _parse_sacct(data)
    return jobs[job_id]


def get_all(*job_ids):
    """ Get job data for multiple jobs. """
    job_ids = tuple(str(j) for j in job_ids)
    data = _sacct_request(*job_ids)
    jobs = _parse_sacct(data)

    for job_id in job_ids:
        if not job_id in jobs:
            msg = 'Error fetching job data for {}'.format(job_id)
            log.critical(msg)

    return jobs


def scancel(*args):
    cmd_args = ('scancel',) + args

    log.debug('Launching: {}'.format(cmd_args))
    return subprocess.call(cmd_args)


def sbatch(*args, stdin_data=None):
    cmd_args = ('sbatch', '--parsable') + args

    if stdin_data and not isinstance(stdin_data, bytes):
        stdin_data = stdin_data.encode()

    log.debug('Launching: {}'.format(cmd_args))
    p = subprocess.Popen(
        cmd_args,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE
    )

    stdout, stderr = p.communicate(input=stdin_data)

    if p.returncode != 0:
        raise ChildProcessError(str(vars(p)))

    jid, _, cluster = stdout.decode().strip().partition(';')
    log.critical("Submitted batch job {}".format(jid))
    return Sopen(jid, cluster=cluster)

