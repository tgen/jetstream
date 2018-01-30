""" Async runners operate on workflows to execute tasks. They make use
of workflow methods .__next__() and .__send__() """

import logging
import shlex
import asyncio, sys

from jetstream import plugins

log = logging.getLogger(__name__)

if sys.platform == "win32":
    # On Windows, the default event loop is SelectorEventLoop
    # which does not support subprocesses.
    loop = asyncio.ProactorEventLoop()
    asyncio.set_event_loop(loop)
else:
    loop = asyncio.get_event_loop()


async def spawn(wf, task, **kwargs):
    """ Spawn a subprocess, this returns a coroutine """

    # Capture stdout/stderr by default. This may cause memory issues 
    # if the subprocesses produce too much stream data. but less code for now.
    try:
        if 'stdout' not in kwargs:
            kwargs['stdout'] = asyncio.subprocess.PIPE
        if 'stderr' not in kwargs:
            kwargs['stderr'] = asyncio.subprocess.PIPE

        node_name, node_data = task

        # Prep command for execution on slurm
        cmd = slurm(node_name)
        cmd_args = shlex.split(cmd)

        proc = await asyncio.create_subprocess_exec(*cmd_args, **kwargs)
        log.critical('Spawned: {} pid: {}'.format(cmd, proc.pid))

        stdout, stderr = await proc.communicate()

        record = {
            "pid": proc.pid,
            "returncode": proc.returncode,
            "stdout": stdout,
            "stderr": stderr,
        }

        wf.__send__((node_name, record))

    except Exception as err:
        loop.stop()
        raise err


async def wf_walker(wf):
    """ Walks a workflow, submitting tasks as they become ready. 
    This returns a coroutine. """
    while 1:
        try:
            task = next(wf)
        except StopIteration:
            log.critical('wf raised stop iteration')
            break

        if task is None:
            await asyncio.sleep(1)
        else:
            loop.create_task(spawn(wf, task))


def run(wf, debug=True):
    log.critical('starting event loop')
    loop.set_debug(debug)
    loop.run_until_complete(wf_walker(wf))
    loop.close()
    log.critical('wf appears to be done')


async def _spawn(plugin_id, cb):
    try:
        stdout = stderr = asyncio.subprocess.PIPE

        cmd = plugins[plugin_id]
        cmd_args = shlex.split(cmd)

        proc = await asyncio.create_subprocess_exec(*cmd_args, stdout=stdout, stderr=stderr)
        log.critical('Spawned: {} pid: {}'.format(cmd, proc.pid))

        stdout, stderr = await proc.communicate()

        record = {
            "pid": proc.pid,
            "returncode": proc.returncode,
            "stdout": stdout,
            "stderr": stderr,
        }

        cb((plugin_id, record))

    except Exception as err:
        loop.stop()
        raise err


def run_one(plugin_id, debug=True):
    log.critical('starting event loop')
    loop.set_debug(debug)
    loop.run_until_complete(_spawn(plugin_id, print))
    log.critical('event loop stopped')

