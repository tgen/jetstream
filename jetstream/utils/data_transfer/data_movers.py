import paramiko
import time
import itertools
import logging

log = logging.getLogger('jetstream_central')


class Wrangler(object):
    """Manages a set of data movers as a simple round-robin load balancer"""

    def __init__(self, hosts):
        self.hosts = hosts
        self.cycler = itertools.cycle(self.hosts)

    def next(self):
        """Return the next data mover."""
        return next(self.cycler)


class LowestLoadWrangler(Wrangler):
    """Manages a set of data movers either by finding the machine with the
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
    """Manages a set of data movers either by finding the machine with the
    lowest current load"""

    def __init__(self, hosts, max_rsyncs=8, delay=10, timeout=None):
        super(RsyncLimitedWrangler, self).__init__(hosts)
        self.max_rsyncs = max_rsyncs
        self.delay = delay
        self.timeout = timeout

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

            if self.timeout and elapsed >= self.timeout:
                break

            for host in self.hosts:
                if check_rsyncs(host) < self.max_rsyncs:
                    return host

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
    log.info('Checking rsyncs on: {}'.format(host))

    ssh = paramiko.SSHClient()
    ssh.load_system_host_keys()
    ssh.connect(host)

    stdin, stdout, stderr = ssh.exec_command('ps aux | grep rsync')
    _ = stdout.channel.recv_exit_status()
    rsyncs = stdout.read().decode('utf8').splitlines()

    return max(len(rsyncs) - 2, 0)
