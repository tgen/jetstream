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
active_states = {'CONFIGURING', 'COMPLETING', 'RUNNING', 'SPECIAL_EXIT',}
inactive_states = {'BOOT_FAIL', 'CANCELLED', 'COMPLETED', 'FAILED',
                   'NODE_FAIL', 'PENDING', 'PREEMPTED', 'REVOKED',
                   'STOPPED', 'SUSPENDED', 'TIMEOUT',}
failed_states = {'BOOT_FAIL', 'CANCELLED', 'FAILED', 'NODE_FAIL',}
completed_states = {'COMPLETED',}


class SlurmJob(object):
    # TODO this should be a mix-in that can be used by a generic job class
    def __init__(self, props):
        self.properties = props

    def __getattr__(self, item):
        return self.properties[item]

    def __repr__(self):
        return "SlurmJob(%s)" % self.JobID

    def __str__(self):
        return json.dumps(self.properties)

    @property
    def is_active(self):
        if self.State in active_states:
            return True
        else:
            return False

    @property
    def is_failed(self):
        if self.State in failed_states:
            return True
        else:
            return False

    @property
    def is_complete(self):
        if self.State in completed_states:
            return True
        else:
            return False

    def update(self):
        # TODO This class could be useful if exposed with an update feature
        raise NotImplementedError


def job_ids_to_arguments(job_ids):
    """ Given list of job ids, return as list of arguments with each job
    id prefixed with -j """
    job_ids_prefixed = ' '.join(['-j %s' % jid for jid in job_ids])
    job_ids_split = shlex.split(job_ids_prefixed)
    return job_ids_split


def query_sacct(job_ids):
    """ Run sacct query for given job_ids and returns a list of records """
    log.debug('query_sacct: %s' % job_ids)
    if isinstance(job_ids, (str, int)):
        job_ids = (job_ids,)

    # Build the query and execute sacct
    job_id_args = job_ids_to_arguments(job_ids)
    cmd_args = ['sacct', '-XP', '--format', 'all'] + job_id_args

    log.debug('Launching: %s' % ' '.join(cmd_args))
    res = subprocess.check_output(cmd_args)

    # With python3.x res will be a bytes object that we need to decode
    try:
        res = res.decode()
    except AttributeError:
        pass

    # Convert the sacct report to an object
    records = []
    lines = res.splitlines()
    header = lines.pop(0)
    columns = header.strip().split('|')
    for line in lines:
        row = OrderedDict(zip(columns, line.strip().split('|')))
        records.append(row)
    return records


def get_jobs(job_ids):
    """ Run sacct query for given job_ids, returns list of SlurmJobs """
    records = query_sacct(job_ids)
    jobs = []
    for record in records:
        jobs.append(SlurmJob(record))
    return jobs
