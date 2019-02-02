"""Copy projects to/from isilon storage partitions"""
import os
import sys
import paramiko
import time
import itertools
import argparse
import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urlparse

log = logging.getLogger(__name__)
IN_PROGRESS_FILENAME = 'transfer_in_progress...'


class Wrangler(object):
    """Manages a set of data movers as a simple round-robin load balancer"""
    def __init__(self, hosts):
        self.hosts = hosts
        self.cycler = itertools.cycle(self.hosts)

    def next(self):
        """Return the next data mover."""
        return next(self.cycler)


class LowestLoadWrangler(Wrangler):
    """Manages a set of data movers by finding the machine with the
    lowest current load"""
    def next(self):
        """This will search all hosts for the one with the lowest current load.
        SSH Errors will fall back to round-robin selection of host."""
        try:
            return self.lowest_load()
        except paramiko.SSHException as e:
            log.exception(e)

        return next(self.cycler)

    def lowest_load(self):
        """Finds the host with the lowest load"""
        load = {h: check_load(h) for h in self.hosts}
        lowest = min(load, key=load.get)
        return lowest


class RsyncLimitedWrangler(Wrangler):
    """Manages a set of data movers by finding the machine with less than the
    set number of rsync processes. This will delay until a data mover is found
    up to a set timeout, or forever if timeout is None"""
    def __init__(self, hosts, max_rsyncs=3, delay=10, timeout=None):
        super(RsyncLimitedWrangler, self).__init__(hosts)
        self.max_rsyncs = max(int(max_rsyncs), 0)
        self.delay = max(int(delay), 0)

        if timeout is None:
            self.timeout = timeout
        else:
            self.timeout = max(int(timeout), 0)

    def next(self):
        """This will search all hosts for the one with the lowest current load.
        SSH Errors will fall back to round-robin selection of host."""
        try:
            return self.find_open_dm()
        except paramiko.SSHException as e:
            log.exception(e)

        return next(self.cycler)

    def find_open_dm(self):
        """Finds the host with the lowest load"""
        start = time.time()

        while 1:
            elapsed = time.time() - start

            for host in self.cycler:
                if check_rsyncs(host) < self.max_rsyncs:
                    return host

            if self.timeout is not None and elapsed >= self.timeout:
                break
            else:
                time.sleep(self.delay)

        log.warning('Search timed out!')
        return next(self.cycler)


def check_load(host, period=2):
    """Check load on a host machine.
    Period can be any integer between 0 and 2 corresponding to time intervals:
    last 1 minute, last 5 minutes, or last 10 minutes """
    log.info('Checking load on: {}'.format(host))

    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.connect(host)

    stdin, stdout, stderr = ssh.exec_command('cat /proc/loadavg')
    _ = stdout.channel.recv_exit_status()
    loadavg = stdout.read().decode('utf8').split()

    ssh.close()

    load = loadavg[int(max(0, min(period, 2)))]
    return float(load)


def check_rsyncs(host):
    """Check for rsync jobs on a host machine."""
    log.info('Checking rsyncs on: {}'.format(host))

    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.connect(host)

    stdin, stdout, stderr = ssh.exec_command('ps aux | grep rsync')
    _ = stdout.channel.recv_exit_status()
    rsync_lines = stdout.read().decode('utf8').splitlines()
    n_rsyncs = max(len(rsync_lines) - 2, 0)  # Adjust for the current grep and ssh

    ssh.close()

    log.info(f'{n_rsyncs} rsyncs on {host}!')
    return n_rsyncs


def rsync(src_path, dest_path, *args):
    cmd = ['rsync', '-e', 'ssh -v -o "BatchMode yes"', ]
    cmd += list(args)
    cmd += [src_path, dest_path]
    log.critical('Launching: {}'.format(' '.join(cmd)))
    return subprocess.run(cmd, check=True)


def rsync_project_to_isilon(src_list, dest, wrangler, *args):
    """ This function will take URI's but only utilizes the path. It assumes
    that the transfer is occuring: host ---> isilon and sets the destination
    hostname based on current data mover traffic. Using this function with
    a uri that is not on isilon will cause problems."""
    dest = urlparse(dest).path  # If its a URI, just give me the path

    for src in src_list:
        if not os.path.isdir(src):
            raise ValueError('src is not a directory')

        if src.endswith('/'):
            src = src.rstrip('/')

        target = os.path.join(dest, os.path.basename(src))

        host = wrangler.next()
        ssh = paramiko.SSHClient()
        ssh.load_system_host_keys()
        ssh.connect(host)

        log.info('Making progress file...')
        cmd = f'mkdir -p $(dirname {target}) && touch {target}/{IN_PROGRESS_FILENAME}'
        stdin, stdout, stderr = ssh.exec_command(cmd)
        stdout.channel.recv_exit_status()

        try:
            rsync(src, '{}:{}'.format(host, dest), '-avu', *args)
        finally:
            log.info('Removing progress file...')
            stdin, stdout, stderr = ssh.exec_command(f'rm {target}/{IN_PROGRESS_FILENAME}')
            stdout.channel.recv_exit_status()
            ssh.close()


def rsync_project_from_isilon(src_list, dest, wrangler, *args):
    """ This function will take URI's but only utilizes the path. It assumes
    that the transfer is occuring: isilon ---> host and sets the destination
    hostname based on current data mover traffic. Using this function with
    a uri that is not on isilon will cause problems. """
    if not os.path.isdir(dest):
        raise ValueError('dest is not a directory')

    for src in src_list:
        src = urlparse(src).path  # If its a URI, just give me the path

        if src.endswith('/'):
            src = src.rstrip('/')

        target = os.path.join(dest, os.path.basename(src))
        host = wrangler.next()

        log.info('Making progress file...')
        target_dir = os.path.dirname(target)
        os.makedirs(target_dir, exist_ok=True)
        with open(f'{target_dir}/{IN_PROGRESS_FILENAME}', 'w'):
            pass

        try:
            rsync('{}:{}'.format(host, src), dest, '-avu', *args)
        finally:
            log.info('Removing progress file...')
            try:
                os.remove(f'{target_dir}/{IN_PROGRESS_FILENAME}')
            except FileNotFoundError:
                pass


def rsync_bulk_with_uris(src_list, dest, wrangler, max_workers=3, *args):
    """Parallel rsync with URI support
    This supports netloc for anywhere, whereas t """
    executor = ThreadPoolExecutor(max_workers=max_workers)

    with executor:
        jobs = {}

        for src in src_list:
            src_path = resolve_location(src, wrangler)
            dest_path = resolve_location(dest, wrangler)
            jobs[executor.submit(rsync, src_path, dest_path, *args)] = src

        for future in as_completed(jobs):
            future.result()


def resolve_location(value, wrangler):
    """Given a file path, convert to valid uri, and return a formatted location
    for rsync."""
    if isinstance(value, str):
        uri = urlparse(value, scheme='file')
    else:
        uri = value

    if uri.scheme == 'file' and uri.netloc == 'isilon.tgen.org':
        new_host = wrangler.next()
        uri = uri._replace(netloc=new_host)

    if uri.netloc:
        return '{}:{}'.format(uri.netloc, uri.path)
    else:
        return uri.path


def arg_parser():
    parser = argparse.ArgumentParser(
        description='Transfer project folders to/from Isilon storage partitions.'
                    'Any additional arguments will be passed to rsync.'
    )

    parser.add_argument(
        'action',
        choices=['to_isilon', 'from_isilon', 'bulk']
    )

    parser.add_argument(
        'paths',
        nargs='+'
    )

    parser.add_argument(
        '--data-movers',
        help='Modify the set of data movers to use, This should be a semicolon-'
             'delimited list of host names',
        default='dback-data1.tgen.org;dback-data2.tgen.org;dback-data3.tgen.org'
    )

    parser.add_argument(
        '--max-rsyncs',
        help='Maximum number of rsyncs allowed per data mover.',
        type=int,
        default=3
    )

    parser.add_argument(
        '--max-workers',
        type=int,
        default=1
    )

    return parser


def main(args=None):
    logging.basicConfig(level=logging.INFO)

    parser = arg_parser()
    args, remainder = parser.parse_known_args(args)
    log.info(f'Cmd args: {args}\nRsync args: {remainder}')

    if len(args.paths) < 2:
       raise ValueError('No destination given')
    else:
        src_list = args.paths[:-1]
        dest = args.paths[-1]

    wrangler = RsyncLimitedWrangler(
        hosts=args.data_movers.split(';'),
        max_rsyncs=args.max_rsyncs
    )

    if args.action == 'to_isilon':
        rsync_project_to_isilon(src_list, dest, wrangler, *remainder)
    elif args.action == 'from_isilon':
        rsync_project_from_isilon(src_list, dest, wrangler, *remainder)
    elif args.action == 'bulk':
        rsync_bulk_with_uris(src_list, dest, wrangler, args.max_workers, *remainder)


def to_isilon():
    """Entry point helper"""
    sys.argv.insert(1, 'to_isilon')
    main()


def from_isilon():
    """Entry point helper"""
    sys.argv.insert(1, 'from_isilon')
    main()


def bulk():
    """Entry point helper"""
    sys.argv.insert(1, 'bulk')
    main()


if __name__ == '__main__':
    main()

