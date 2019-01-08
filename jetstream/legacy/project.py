import fnmatch
import logging
import re
from os import path, walk, listdir

from jetstream.backends import slurm

log = logging.getLogger(__name__)


class Project(object):
    def __init__(self, project):
        self.path = path.realpath(project)
        self.name = path.basename(self.path)
        self.short_name = re.subn('_ps\d*$', '', self.name)[0]
        #self.start_date = re.match('.*_ps(\d*)').groups()[0]

        # Check for a config file
        self.config_path = path.join(self.path, self.short_name + '.config')
        if not path.exists(self.config_path):
            raise FileNotFoundError('Config file not found: {}'.format(
                self.config_path))

    def __repr__(self):
        return "Project({})".format(self.name)

    @property
    def is_complete(self):
        """ Faster than check_status """
        return path.exists(path.join(self.path, 'project.finished'))

    def check_status(self):
        """ There is no single source of truth when it comes to a project
        status. The best we can do is make a guess based on:

        1) Whether or not a "project.finished" file exists
        2) .Failed files exist in the project directory
        3) .Queue files exist in the project directory
        4) Checking the state of job ids recorded in the project logs

        This method only takes into account 1 and 2. Combined this
        with Project.get_jobs() and you can get a full picture on the
        status of a project """
        if path.exists(path.join(self.path, 'project.finished')):
            return 'complete'
        elif find_failed_signals(self.path):   # This operation can take a while
            return 'failed'
        else:
            return 'incomplete'

    def get_jobs(self):
        # Using this pattern in order to batch request job info from
        # sacct. We could just iterate over jobs calling job.update(),
        # but each call to job.update() starts a separate sacct process.
        # It's faster to batch them together in a single request
        jids = list(self._jids())
        if jids:
            return slurm.sacct(*jids)
        else:
            return []

    def _logs(self):
        log_dir = path.join(self.path, 'logs/')
        for l in listdir(log_dir):
            log_path = path.join(log_dir, l)
            if path.isfile(log_path) and log_path.endswith('.txt'):
                yield log_path
        raise StopIteration

    def _jids(self):
        for l in self._logs():
            for jid in find_jids_in_log(l):
                yield jid
        raise StopIteration

    def report(self, fast=False, all_jobs=False, col_size=20):
        if fast:
            if self.is_complete:
                status = 'complete'
            else:
                status = 'incomplete'

            rep = "%s - %s\n" % (self.name, status)
            rep += "%s\n" % self.path

        else:
            rep = "%s - %s\n" % (self.name, self.check_status())
            rep += "%s\n" % self.path

            # Build the header
            h_id = "ID".ljust(col_size)
            h_name = "Name".ljust(col_size)
            h_state = "State".ljust(col_size)
            h_start = "Start".ljust(col_size)
            h_end = "End".ljust(col_size)
            h_elapsed = "Elapsed".ljust(col_size)
            header_values = (h_id, h_name, h_state, h_start, h_end, h_elapsed)
            rep += " ".join(header_values) + '\n'


            for j in self.get_jobs():
                if not all_jobs and j.is_done() and j.is_ok():
                    continue
                
                d = j.job_data  # Use the sacct data for each job
                if d is None:
                    print(f'Error loading job data for: {j}')
                    continue
                try:
                    job_id = d['JobID'][:12].ljust(12)
                    job_name = d['JobName'][:col_size].ljust(col_size)
                    job_state = d['State'][:col_size].ljust(col_size)
                    job_start = d['Start'][:col_size].ljust(col_size)
                    job_end = d['End'][:col_size].ljust(col_size)
                    job_elapsed = d['Elapsed'][:col_size].ljust(col_size)
                    values = (job_id, job_name, job_state, job_start, job_end,
                              job_elapsed)

                    rep += " ".join(values) + '\n'
                except KeyError:
                    print(f'Error loading job data for: {j}')
                    continue

        return rep


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


def find_jids_in_log(log_path):
    """ Given path to a log file, search the log file for job ids """
    log.debug('Searching for job ids in: %s' % log_path)

    pat = re.compile("Submitted batch job (\d*)")

    with open(log_path, 'r') as fp:
        for line in fp.readlines():
            match = pat.match(line)
            if match:
                yield match.groups()[0]

    raise StopIteration
