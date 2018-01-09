""" Async runners operate on workflows to execute tasks. They make use
of workflow methods .__next__() and .__send__() """

import asyncio, sys
from jetstream.analysis import launcher


if sys.platform == "win32":
    # On Windows, the default event loop is SelectorEventLoop
    # which does not support subprocesses.
    loop = asyncio.ProactorEventLoop()
    asyncio.set_event_loop(loop)
else:
    loop = asyncio.get_event_loop()


async def spawn(wf, task, *args, **kwargs):
    """ Spawn a subprocess, returns a coroutine """

    # Capture stdout/stderr by default. This may cause memory issues 
    # if the subprocesses produce too much stream data. but less code for now.
    if 'stdout' not in kwargs:
        kwargs['stdout'] = asyncio.subprocess.PIPE
    if 'stderr' not in kwargs:
        kwargs['stderr'] = asyncio.subprocess.PIPE

    node_name, node_data = task
    cmd = launcher.faux_cmd(node_name)
    proc = await asyncio.create_subprocess_exec(*cmd, *args, **kwargs)
    print('Spawned:', str(cmd), ' pid:', str(proc.pid))
    stdout, stderr = await proc.communicate()
    wf.__send__((node_name, proc))


async def wf_walker(loop, wf):
    """ Walks a workflow, submitting tasks as they become ready. 
    This returns a coroutine. """
    while 1:
        try:
            task = next(wf)
        except StopIteration:
            print('wf raised stop iteration')
            break

        if task is None:
            await asyncio.sleep(1)
        else:
            loop.create_task(spawn(wf, task))


def run(wf, debug=True):
    print('starting event loop')
    loop.set_debug(debug)
    loop.run_until_complete(wf_walker(loop, wf))
    loop.close()
    print('wf appears to be done')

