import re
import fnmatch
import logging
from os import path, walk, listdir

from jetstream.batch_schedulers.slurm import submission_pattern, get_jobs

log = logging.getLogger(__name__)


class Project(object):
    def __init__(self, project, deep=False):
        # Resolve project path information
        self.path = path.realpath(project)
        self.name = path.basename(self.path)
        self.short_name = re.subn('_ps\d*$', '', self.name)[0]

        # Load the config file
        self.config_path = path.join(self.path, self.short_name + '.config')
        if not path.exists(self.config_path):
            raise AttributeError('Config file no found: %s' % self.config_path)

        # Check the status
        self.fail_files = None
        self.queue_files = None
        self.status = 'incomplete'
        self.update(deep=deep)

    def __repr__(self):
        return "Project(%s, %s)" % (self.name, self.status)

    def serialize(self):
        """ Returns a dictionary with the project attributes and all jobs """
        attrs = {k: v for k, v in self.__dict__.items() if
                 not k.startswith('_')}
        attrs['jobs'] = [j.serialize() for j in self.jobs()]
        return attrs

    def update(self, deep=False):
        """ There is no single source of truth when it comes to a project's
        status. The best we can do is make a guess based on:

        1) Whether or not a "project.finished" file exists
        2) .Failed or .Queue files exist in the project directory
        3) Checking the state of job ids recorded in the project logs

        This method attempts to turn all of those clues into a single status
        value """
        if path.exists(path.join(self.path, 'project.finished')):
            self.status = 'complete'

        if deep:
            # These operations can take quite a while, so they're optional
            self.fail_files = find_failed_signals(self.path)
            self.queue_files = find_queued_signals(self.path)

        if self.fail_files:
            self.status = 'failed'
        elif self.queue_files:
            self.status = 'active'

        self.update_jobs()

        return self.status

    def update_jobs(self):
        # Using this pattern here so that we batch request job info from
        # sacct. We could just iterate over jobs calling job.update(),
        # but each call to job.update() starts an sacct process. It's faster
        # to batch them
        jids = self._jids()
        self.jobs = get_jobs(*jids)

    def _logs(self):
        log_dir = path.join(self.path, 'logs/')
        for log in listdir(log_dir):
            log_path = path.join(log_dir, log)
            if path.isfile(log_path) and log_path.endswith('.txt'):
                yield log_path
        raise StopIteration

    def _jids(self):
        for log in self._logs():
            for jid in find_jids_in_log(log):
                yield jid
        raise StopIteration

    @property
    def active_jobs(self):
        active_jobs = [j for j in self.jobs if j.is_active]
        return active_jobs

    @property
    def complete_jobs(self):
        active_jobs = [j for j in self.jobs if j.is_complete]
        return active_jobs

    @property
    def failed_jobs(self):
        active_jobs = [j for j in self.jobs if j.is_failed]
        return active_jobs

    @property
    def is_complete(self):
        if self.status == 'complete':
            return True
        else:
            return False

    @property
    def is_active(self):
        if self.status != 'complete':
            return True
        else:
            return False


def find_failed_signals(project_dir):
    fails = []
    for root, dirnames, filenames in walk(project_dir):
        for filename in fnmatch.filter(filenames, '*Fail'):
            fails.append(path.join(root, filename))
    return fails


def find_queued_signals(project_dir):
    queues = []
    for root, dirnames, filenames in walk(project_dir):
        for filename in fnmatch.filter(filenames, '*Queue'):
            queues.append(path.join(root, filename))
    return queues


def find_jids_in_log(log_path, pat=submission_pattern):
    """ Given path to a log file, search the log file for job ids """
    log.debug('Searching for job ids in: %s' % log_path)

    with open(log_path, 'r') as fp:
        for line in fp.readlines():
            match = pat.match(line)
            if match:
                yield match.groups()[0]

    raise StopIteration


def build_plain_text_report(project, col_size=20):
    rep = "%s - %s\n" % (project.name, project.status)
    rep += "%s\n" % project.path

    # Build the header
    h_id = "ID".ljust(col_size)
    h_name = "Name".ljust(col_size)
    h_state = "State".ljust(col_size)
    h_start = "Start".ljust(col_size)
    h_end = "End".ljust(col_size)
    h_elapsed = "Elapsed".ljust(col_size)
    header_values = (h_id, h_name, h_state, h_start, h_end, h_elapsed)
    rep += " ".join(header_values) + '\n'

    for j in project.jobs():
        d = j.sacct  # Use the sacct data for each job
        job_id = d['JobID'][:12].ljust(12)
        job_name = d['JobName'][:col_size].ljust(col_size)
        job_state = d['State'][:col_size].ljust(col_size)
        job_start = d['Start'][:col_size].ljust(col_size)
        job_end = d['End'][:col_size].ljust(col_size)
        job_elapsed = d['Elapsed'][:col_size].ljust(col_size)
        values = (job_id, job_name, job_state, job_start, job_end, job_elapsed)
        rep += " ".join(values) + '\n'
    return rep
