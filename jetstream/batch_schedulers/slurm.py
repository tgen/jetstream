import re
import json
import shlex
import subprocess
import logging
from collections import OrderedDict

log = logging.getLogger(__name__)

submission_pattern = re.compile("Submitted batch job (\d*)")
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
                 'PENDING',}
inactive_states = {'BOOT_FAIL', 'CANCELLED', 'COMPLETED', 'FAILED',
                   'NODE_FAIL', 'PREEMPTED', 'REVOKED',
                   'STOPPED', 'SUSPENDED', 'TIMEOUT',}
failed_states = {'BOOT_FAIL', 'CANCELLED', 'FAILED', 'NODE_FAIL',}
completed_states = {'COMPLETED',}


class SlurmJob(object):
    def __init__(self, jid, sacct=None, cluster=None):
        self.jid = jid
        self.sacct = sacct
        self.cluster = cluster

    def __repr__(self):
        return "SlurmJob(%s)" % self.jid

    def __str__(self):
        return json.dumps(self.__dict__)

    def serialize(self):
        return self.__dict__

    @property
    def is_active(self):
        if self.sacct['State'] in active_states:
            return True
        else:
            return False

    @property
    def is_failed(self):
        if self.sacct['State'] in failed_states:
            return True
        else:
            return False

    @property
    def is_complete(self):
        if self.sacct['State'] in completed_states:
            return True
        else:
            return False

    def update(self):
        sacct_data = query_sacct(self.jid)
        matches = [r for r in sacct_data if r['JobID'] == self.jid]

        if len(matches) > 1:
            msg = "Sacct returned more than one record for {}".format(self.jid)
            raise ValueError(msg)
        elif len(matches) < 1:
            msg = "Sacct returned less than one record for {}".format(self.jid)
            raise ValueError(msg)
        else:
            self.sacct = matches[0]

        return self.sacct


def query_sacct(*job_ids, all=False):
    """ Run sacct query for given job_ids and returns a list of records """
    log.debug('query_sacct: {}'.format(str(job_ids)))

    cmd_prefix = ['sacct', '-XP', '--format', 'all']

    if job_ids:
        job_ids = ' '.join(['-j %s' % jid for jid in job_ids if jid])
        cmd_args = cmd_prefix + shlex.split(job_ids)
    else:
        if not all:
            raise ValueError('must give job ids or specify all=True')
        cmd_args = cmd_prefix

    log.debug('Launching: %s' % ' '.join(cmd_args))
    res = subprocess.check_output(cmd_args).decode()

    # Convert the sacct report to an object
    records = []
    lines = res.splitlines()
    header = lines.pop(0).strip().split('|')

    for line in lines:
        row = OrderedDict(zip(header, line.strip().split('|')))
        records.append(row)

    return records


def get_jobs(*args, **kwargs):
    """ Run batch query for slurm jobs, returns a list of SlurmJobs """
    jobs = []
    records = query_sacct(*args, **kwargs)
    for rec in records:
        s = SlurmJob(rec['JobID'], sacct=rec)
        jobs.append(s)

    return jobs


def srun(*args):
    """ Srun a command, additional args are added to srun prefix """
    cmd_args = ('srun',) + args
    log.debug('launching: {}'.format(cmd_args))
    return subprocess.check_output(cmd_args)


def sbatch(*args, stdin_data=None):
    cmd_args = ('sbatch', '--parsable') + args

    log.debug('launching: {}'.format(cmd_args))
    p = subprocess.Popen(
        cmd_args,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE
    )

    stdout, stderr = p.communicate(input=stdin_data)

    if p.returncode != 0:
        raise ChildProcessError(str(vars(p)))

    jid, _, cluster = stdout.decode().partition(';')
    return SlurmJob(jid, cluster=cluster)


