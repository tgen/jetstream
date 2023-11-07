import asyncio
import glob
import itertools
import json
import logging
import os
import sys
import re
import shlex
import shutil
import signal
import subprocess
import time
from asyncio import Lock, BoundedSemaphore, create_subprocess_shell, CancelledError
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
SLURM_SOLO_OPTIONS = settings['slurm_solo_options'].get(list)
SLURM_PRESETS = settings['slurm_presets'].get(dict)
RUNNER_PRESETS = settings['singularity_presets'].get(dict)


class SlurmSingularityBackend(BaseBackend):
    """SlurmSingularityBackend will spawn tasks using a Slurm batch scheduler.

    The spawn coroutine will return when the slurm job ID for a task is
    complete. This works by maintaining a dict of SlurmBatchJobs, and
    periodically asking for updates from sacct."""
    count = itertools.count()
    respects = ('cmd', 'stdin', 'stdout', 'stderr', 'cpus', 'mem', 'walltime',
                'slurm_args')

    def __init__(self,
                 sacct_frequency=60,
                 sbatch_delay=0.5,
                 sbatch_executable=None,
                 sbatch_account=None,
                 sacct_fields=('JobID', 'Elapsed', 'State', 'ExitCode'),
                 job_monitor_max_fails=5,
                 max_jobs=1024,
                 singularity_executable=None,
                 input_file_validation=False):
        """SlurmSingularityBackend submits tasks as jobs to a Slurm batch cluster

        :param sacct_frequency: Frequency in seconds that job updates will
        be requested from sacct
        :param sbatch: path to the sbatch binary if not on PATH
        """
        super(SlurmSingularityBackend, self).__init__()
        self.sacct_frequency = sacct_frequency
        self.sbatch_delay = sbatch_delay
        self.sbatch_executable = sbatch_executable
        self.sbatch_account = sbatch_account
        self.sacct_fields = sacct_fields
        self.sbatch_lock = Lock()
        self.job_monitor_max_fails = job_monitor_max_fails
        self.max_jobs = max_jobs
        self.jobs = dict()
        self.input_file_validation = input_file_validation

        self.coroutines = (self.job_monitor,)
        self._next_update = datetime.now()

        if self.sbatch_executable is None:
            self.sbatch_executable = shutil.which('sbatch') or 'sbatch'

        # with open(os.devnull, 'w') as devnull:
        #     subprocess.run(
        #         [self.sbatch_executable, '--version'],
        #         check=True,
        #         stdout=devnull,
        #         stderr=devnull
        #     )

        self.singularity_executable = singularity_executable
        if self.singularity_executable is None:
            self.singularity_executable = shutil.which('singularity') or 'singularity'

        # To ensure pulls have exclusive use of singularity
        self._singularity_run_sem = BoundedSemaphore(self.max_jobs)
        self._singularity_pull_lock = Lock()
        self._singularity_pull_cache = {}

        signal.signal(signal.SIGABRT, self.cancel)
        signal.signal(signal.SIGHUP, self.cancel)
        signal.signal(signal.SIGIOT, self.cancel)
        signal.signal(signal.SIGQUIT, self.cancel)
        signal.signal(signal.SIGTERM, self.cancel)
        signal.signal(signal.SIGINT, self.cancel)

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

        input_filenames = task.directives.get('input', [])
        output_filenames = task.directives.get('output', [])

        container = task.directives.get('container', None)
        digest = task.directives.get('digest', None)
        if container is None:
            # raise RuntimeError(f'container argument missing for task: {task.name}')
            log.debug(f'Task: {task.name} is missing a container definition! Using basic slurm submission')
            singularity_image = None
            singularity_hostname = None
            docker_authentication_token = None
        else:
            if os.path.exists(container):
                singularity_image = container
                singularity_hostname = None
                docker_authentication_token = None
            else:
                try:
                    image, tag = container.split(':')
                except ValueError:
                    log.debug(f'Tag not defined for {container}, assuming latest')
                    image = container
                    tag = 'latest'

                if digest is None:
                    singularity_image = f"docker://{image}:{tag}"
                else:
                    # Stripping sha256 in case it was already included in digest
                    digest = re.sub('^sha256:', '', digest)
                    singularity_image = f"docker://{image}@sha256:{digest}"

                singularity_hostname = task.directives.get('singularity_hostname', None)

                docker_authentication_token = task.directives.get('docker_authentication_token', None)

                log.debug(f'Task: {task.name}, going to pull: {singularity_image}')
                try:
                    if singularity_image in self._singularity_pull_cache:
                        pass
                    else:
                        async with self._singularity_pull_lock:
                            for i in range(self.max_jobs):
                                await self._singularity_run_sem.acquire()
                            pull_command_run_string = ""
                            if docker_authentication_token is not None:
                                pull_command_run_string += f"""SINGULARITY_DOCKER_USERNAME='$oauthtoken' SINGULARITY_DOCKER_PASSWORD={docker_authentication_token} """
                            pull_command_run_string += f'singularity exec --cleanenv --nohttps {singularity_image} true'
                            log.debug(f'Task: {task.name}, pulling: {pull_command_run_string}')
                            _p = await create_subprocess_shell(pull_command_run_string,
                                                               stdout=asyncio.subprocess.PIPE,
                                                               stderr=asyncio.subprocess.PIPE)
                            _stdout, _stderr = await _p.communicate()
                            log.debug(f'Task: {task.name}, pulled, stdout: {_stdout}')
                            log.debug(f'Task: {task.name}, pulled, stderr: {_stderr}')
                            self._singularity_pull_cache[singularity_image] = singularity_image
                            for i in range(self.max_jobs):
                                self._singularity_run_sem.release()
                except Exception as e:
                    log.warning(f'Exception during singularity prepull: {e}')
                    p = await create_subprocess_shell("exit 1;")

                log.debug(f'Task: {task.name}, pull complete: {singularity_image}')

        async with self.sbatch_lock:
            time.sleep(self.sbatch_delay)
            job = await sbatch(
                cmd=task.directives['cmd'],
                identity=task.identity,
                singularity_image=singularity_image,
                singularity_image_digest=digest,
                singularity_executable=self.singularity_executable,
                singularity_run_sem=self._singularity_run_sem,
                singularity_hostname=singularity_hostname,
                runner_args=task.directives.get('runner_args'),
                runner_preset=task.directives.get('runner_preset'),
                docker_authentication_token=docker_authentication_token,
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
                queue_args=task.directives.get('queue_args'),
                queue_preset=task.directives.get('queue_preset'),
                sbatch_executable=self.sbatch_executable,
                sbatch_account=self.sbatch_account,
                input_file_validation=self.input_file_validation
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

        if sys.version_info > (3, 8):
            job.event = asyncio.Event()
        else:
            job.event = asyncio.Event(loop=self.runner.loop)
        self.jobs[job.jid] = job

        await job.event.wait()
        log.debug(f'{task.name}: job info was updated')

        if self.sacct_fields:
            job_info = {k: v for k, v in job.job_data.items() if k in self.sacct_fields}
            log.debug(f'{task.name} job.job_data: {job.job_data}')
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

        if self.jid not in data:
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
        if job.jid not in data:
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


async def sbatch(cmd, identity, singularity_image, singularity_executable="singularity",
                 singularity_run_sem=None, singularity_hostname=None, singularity_image_digest=None,
                 runner_args=None, runner_preset=None, docker_authentication_token=None, name=None,
                 input_filenames=[], output_filenames=[], stdin=None, stdout=None, stderr=None,
                 tasks=None, cpus_per_task=1, mem="2G", walltime="1h", comment=None,
                 queue_args=None, queue_preset=None, sbatch_executable=None,
                 sbatch_account=None, input_file_validation=False):

    # determine input/output mounts needed
    singularity_mounts = set()
    if input_file_validation:
        for input_filename_glob_pattern in input_filenames:
            input_filenames_glob = glob.glob(input_filename_glob_pattern)
            if len(input_filenames_glob) == 0:
                raise RuntimeError(f'Task {name}: input file(s) do not exist: {input_filename_glob_pattern}')
            for input_filename in input_filenames_glob:
                input_filename = os.path.abspath(input_filename)
                input_filename_head, input_filename_tail = os.path.split(input_filename)
                if not input_filename_head.startswith(os.getcwd()):
                    singularity_mounts.add(input_filename_head)
    else:
        for input_filename in input_filenames:
            input_filename = os.path.abspath(input_filename)
            input_filename_head, input_filename_tail = os.path.split(input_filename)
            if not input_filename_head.startswith(os.getcwd()):
                singularity_mounts.add(input_filename_head)
    for output_filename in output_filenames:
        output_filename = os.path.abspath(output_filename)
        output_filename_head, output_filename_tail = os.path.split(output_filename)
        if not output_filename_head.startswith(os.getcwd()):
            singularity_mounts.add(output_filename_head)
        os.makedirs(output_filename_head, exist_ok=True)

    mount_strings = []
    for singularity_mount in singularity_mounts:
        mount_strings.append("--bind %s" % (singularity_mount))
    singularity_mounts_string = " ".join(mount_strings).strip()

    # create cmd script
    os.makedirs("jetstream/cmd", mode=0o777, exist_ok=True)
    millis = int(round(time.time() * 1000))
    if name is None:
        name = "script"
    cmd_script_filename = f"jetstream/cmd/{millis}_{identity}.cmd"
    cmd_script_filename = os.path.abspath(cmd_script_filename)
    with open(cmd_script_filename, "w") as cmd_script:
        cmd_script.write(cmd)

    sbatch_args = []

    if name:
        sbatch_args.extend(['-J', name])

    if stdin:
        sbatch_args.extend(['--input', stdin])

    if stdout:
        sbatch_args.extend(['-o', stdout])

    if stderr:
        sbatch_args.extend(['-e', stderr])

    if tasks:
        sbatch_args.extend(['-n', tasks])

    if cpus_per_task:
        sbatch_args.extend(['-c', cpus_per_task])

    if mem:
        sbatch_args.extend(['--mem', mem])
    elif cpus_per_task:
        sbatch_args.extend(['--mem', f"{cpus_per_task*2}G"])
    else:
        sbatch_args.extend(['--mem', "2G"])

    if walltime:
        sbatch_args.extend(['-t', walltime])

    if sbatch_account:
        sbatch_args.extend(['--account', sbatch_account])

    if comment:
        fixed_comment = '"' + comment.replace('"', "'") + '"'
        sbatch_args.extend(['--comment', fixed_comment])

    if queue_args:
        if isinstance(queue_args, str):
            sbatch_args.append(queue_args)
        else:
            sbatch_args.extend(queue_args)

    """
    SLURM_PRESETS:
        The planned usage here is to configure the jetstream settings with a set of presets. This is useful
        in the case of running a pipeline in two different slurm environments. For example:
        Slurm env 1:
            defq = CPU nodes with 7 day walltime
            overflow = CPU nodes with 1 day walltime
        Slurm env 2:
            defq = CPU nodes with 7 day walltime
            scavenge = CPU nodes with 1 day walltime

    Instead of pushing out different implementations of the specific pipeline, we can define the changes
    pseudo-globally. In other words they are slurm environment specific but visually we can maintain the
    same codebase for the pipelines.

    Example implementation:
        "data_mover" = ['-p', 'overflow,data-mover']
        "gpu" = ['-p', 'gpu', '--gres', 'gpu:1']
    """
    if queue_preset:
        sbatch_args.extend(SLURM_PRESETS[queue_preset])

    sbatch_script = "#!/bin/bash\n"
    skip_next = False
    for i in range(0, len(sbatch_args)):
        if skip_next:
            skip_next = False
            continue
        elif sbatch_args[i] in SLURM_SOLO_OPTIONS:
            sbatch_script += f"#SBATCH {sbatch_args[i]}\n"
        elif "=" in sbatch_args[i]:
            sbatch_script += f"#SBATCH {sbatch_args[i]}\n"
        else:
            sbatch_script += f"#SBATCH {sbatch_args[i]} {sbatch_args[i+1]}\n"
            skip_next = True

    singularity_args = []

    """
    RUNNER_PRESETS:
        Similar to slurm preset above, but instead these are passed to singularity
    """
    if runner_preset:
        singularity_args.extend(RUNNER_PRESETS[runner_preset])

    if runner_args:
        if isinstance(runner_args, str):
            singularity_args.append(runner_args)
        else:
            singularity_args.extend(runner_args)

    singularity_exec_args = "--bind $PWD --pwd $PWD --workdir /tmp --cleanenv --contain"
    if os.getenv('JS_PIPELINE_PATH') is not None:
        singularity_exec_args += " --bind {}".format(os.getenv('JS_PIPELINE_PATH'))
        sbatch_script += "export SINGULARITYENV_JS_PIPELINE_PATH={}\n".format(os.getenv('JS_PIPELINE_PATH'))

    if any('gpu' in s for s in [singularity_args, sbatch_args]):
        if all('--nv' not in s for s in singularity_args):
            singularity_exec_args += ' --nv'

    for arg in singularity_args:
        singularity_exec_args += f" {arg}"

    singularity_hostname_arg = ""
    if singularity_hostname is not None:
        singularity_hostname_arg = f"--hostname {singularity_hostname} "

    singularity_run_env_vars = ""
    if docker_authentication_token is not None:
        singularity_run_env_vars += f"""SINGULARITY_DOCKER_USERNAME='$oauthtoken' SINGULARITY_DOCKER_PASSWORD={docker_authentication_token} """

    if singularity_image:
        # CUDA_VISIBLE_DEVICES is a standard method for declaring which GPUs a user is authorized to use - recognized by tensorflow for example
        sbatch_script += f"[[ -v CUDA_VISIBLE_DEVICES ]] && export SINGULARITYENV_CUDA_VISIBLE_DEVICES=\"$CUDA_VISIBLE_DEVICES\"\n"
        # We set the SINGULARITY_CACHEDIR to the default if it isn't defined by the user
        sbatch_script += f"[[ -v SINGULARITY_CACHEDIR ]] || SINGULARITY_CACHEDIR=$HOME/.singularity\n"
        # Searching for the cached image and using it if it exists
        sbatch_script += f"for file in $(find $SINGULARITY_CACHEDIR/cache/oci-tmp -type f); do\n"
        sbatch_script += f"  {singularity_executable} inspect $file 2> /dev/null | grep -q '{singularity_image.split('docker://')[-1]}' && IMAGE_PATH=$file && break\n"
        sbatch_script += f"done\n"
        sbatch_script += f"if [[ -v IMAGE_PATH ]] ; then\n"
        sbatch_script += f"  {singularity_run_env_vars}{singularity_executable} exec {singularity_exec_args} {singularity_hostname_arg}{singularity_mounts_string} $IMAGE_PATH bash {cmd_script_filename}\n"
        sbatch_script += f"else\n"
        sbatch_script += f"  {singularity_run_env_vars}{singularity_executable} exec {singularity_exec_args} {singularity_hostname_arg}{singularity_mounts_string} {singularity_image} bash {cmd_script_filename}\n"
        sbatch_script += f"fi\n"
    else:
        sbatch_script += f"bash {cmd_script_filename}\n"

    if name is None:
        name = "script"
    sbatch_script_filename = f"jetstream/cmd/{millis}_{identity}.sbatch"
    sbatch_script_filename = os.path.abspath(sbatch_script_filename)
    with open(sbatch_script_filename, "w") as sbatch_script_file:
        sbatch_script_file.write(sbatch_script)

    submit_sbatch_args = ["sbatch", sbatch_script_filename]
    remaining_tries = SLURM_SBATCH_RETRY

    while 1:
        if singularity_run_sem is not None:
            await singularity_run_sem.acquire()
        try:
            p = subprocess.run(submit_sbatch_args, stdout=subprocess.PIPE, check=True)
            break
        except subprocess.CalledProcessError:
            if remaining_tries == 0:
                raise
            else:
                remaining_tries -= 1
                log.exception(f'Error during sbatch, retrying in 60s ...')
                time.sleep(60)
        finally:
            if singularity_run_sem is not None:
                singularity_run_sem.release()

    jid = p.stdout.decode().strip().split()[-1]
    job = SlurmBatchJob(jid)
    job.args = submit_sbatch_args
    job.script = sbatch_script_filename
    return job
