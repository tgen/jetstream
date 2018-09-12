import subprocess
import shlex
import paramiko
import logging
from urllib.parse import urlparse

log = logging.getLogger('jetstream_central')


def rsync(src, dest, *args):
    if isinstance(src, str):
        src = urlparse(src)

    if isinstance(dest, str):
        dest = urlparse(dest)

    if src.netloc:
        src_path = '{}:{}'.format(src.netloc, src.path)
    else:
        src_path = src.path

    if dest.netloc:
        dest_path = '{}:{}'.format(dest.netloc, dest.path)
    else:
        dest_path = dest.path

    cmd = ['rsync', '-e', 'ssh -o "BatchMode yes"']
    cmd += list(args) or []
    cmd += [src_path, dest_path]

    log.critical('Launching: {}'.format(cmd))

    return subprocess.call(cmd)


def rsync_via_ssh(src, dest, host, *args):
    """Run an rsync command on a remote host"""
    if isinstance(src, str):
        src = urlparse(src)

    if isinstance(dest, str):
        dest = urlparse(dest)

    if src.netloc:
        src_path = '{}:{}'.format(src.netloc, src.path)
    else:
        src_path = src.path

    if dest.netloc:
        dest_path = '{}:{}'.format(dest.netloc, dest.path)
    else:
        dest_path = dest.path

    cmd = ['rsync']
    cmd += list(args) or []
    cmd += [src_path, dest_path]
    cmd = ' '.join(shlex.quote(arg) for arg in cmd)

    log.critical('Connecting to: {}'.format(host))
    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.connect(host)

    log.critical('Running: {}'.format(cmd))
    stdout = None
    stderr = None

    try:
        stdin, stdout, stderr = ssh.exec_command(cmd, get_pty=True)
        exit_status = stdout.channel.recv_exit_status()
    finally:
        ssh.close()

        if stdout:
            stdout = stdout.read().decode().strip()
            if stdout:
                log.critical(stdout)

        if stderr:
            stderr = stderr.read().decode().strip()
            if stderr:
                log.critical(stderr)

    return exit_status

