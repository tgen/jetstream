import os
import subprocess
from jetstream.tasks import BaseTask


class LocalTask(BaseTask):
    """ Start a task locally with subprocesses.

    BaseTask Directive Handling
    ============================

    cmd
    ----

    Will be launched with subprocess.Popen

    stdin
    ------

    Will be piped to Popen.stdin fd.

    stdout
    -------

    If present, the path will be opened and Popen.stdout will write to the fd.
    Otherwise, will write to: "logs/<task_id>.log"

    stderr
    -------

    If present, the path will be opened and Popen.stderr will write to the fd.
    Otherwise will join stderr with stdout.

    """
    def __init__(self, task_id, task_directives):
        super(LocalTask, self).__init__(task_id, task_directives)
        self._launched = False
        self._proc = None
        self.fds = set()
        self._stdout_fd = None
        self._stderr_fd = None
        self._setup_out_fds()

    def _setup_out_fds(self):
        os.makedirs('logs', exist_ok=True)
        stdout_fd = open(self.stdout_path, 'w')
        self.fds.add(stdout_fd)
        self._stdout_fd = stdout_fd

        if self.stderr_path == self.stdout_path:
            self._stderr_fd = subprocess.STDOUT
        else:
            stderr_fd = open(self.stderr_path, 'w')
            self.fds.add(stderr_fd)
            self._stderr_fd = stderr_fd

    def poll(self):
        return self.proc.poll()

    def kill(self):
        return self.proc.kill()

    def wait(self):
        self.returncode = self.proc.wait()

        for fd in self.fds:
            fd.close()

        return self.returncode

    def launch(self):
        """Launch this task, returns True if launch successful"""
        self._launched = True

        try:
            p = subprocess.Popen(
                self.task_directives.get('cmd') or 'true',
                stdin=subprocess.PIPE,
                stdout=self._stdout_fd,
                stderr=self._stderr_fd,
                shell=True)

            self.extras['pid'] = p.pid
            self.extras['args'] = p.args

            if self.stdin_data is not None:
                if not isinstance(self.stdin_data, bytes):
                    self.stdin_data = self.stdin_data.encode()

                p.stdin.write(self.stdin_data)
                p.stdin.close()

            self._proc = p
            return True

        except BlockingIOError:
            return False
