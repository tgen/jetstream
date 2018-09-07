Developers guide
===================

New backends need to implement a coroutine method "spawn". They can optionally
implement a coroutine "coro" for managing additional background tasks. Here
is an example:

Note that you may want to handle the asyncio.CancelledError in these coroutines
otherwise an error will be reported when the runner exits early. 

.. code_block:: python

    async def anothercoro():
        print('BAHHHHH')


    class Backend(BaseBackend):
        def __init__(self, failure_rate=0.5, max_concurrency=2):
            self.max_concurrency = max_concurrency
            self.failure_rate = failure_rate

        async def coro(self):
            while 1:
                print('AHHHHHHHHHHHH!', self.runner)
                await anothercoro()
                await asyncio.sleep(3)

        async def spawn(self, msg):
            await asyncio.sleep(random.randint(1,3))

            if random.random() < self.failure_rate:
                raise ValueError('uhoh...')

            return 0


Task identity
==============

The Task class defines the objects that are added to the DAG. One goal of this project was 
to build an engine that could identify tasks that have already been run on a project, and
avoid running them again. This was accomplished by assigning an identity to every task, 
which is computed by hashing the task directives. Task state is not included in the identity
so that a "new" task can be compared to a "completed". But, this does create some other 
challenges when combining workflows. 

Problem:

Workflow produces a task that fails:

.. code_block:: python

    >>> t1 = jetstream.Task(cmd='echo hello world && exit 1')
    >>> t1
    <Task b58e7bf6>
  

I run the workflow in my project:

.. code_block:: yaml

    directed: true
    graph:
      jetstream_version: 1.0.0
    links: []
    multigraph: false
    nodes:
    - id: b58e7bf6bc1cfff2437491baea17198d9e43cf7b
      obj:
        cmd: echo hello world && exit 1
        end: '2018-08-29T14:24:25.324076'
        meta: {}
        returncode: 1
        start:
        status: failed
        tid: b58e7bf6bc1cfff2437491baea17198d9e43cf7b


Uh oh, it failed.. I edit the task, run the workflow again:

.. code_block:: python

    >>> t2 = jetstream.Task(cmd='echo hello world && exit 0')
    >>> t2
    <Task c53235a6>
    >>>


My new project workflow:

.. code_block:: yaml

    directed: true
    graph:
      jetstream_version: 1.0.0
    links: []
    multigraph: false
    nodes:
    - id: b58e7bf6bc1cfff2437491baea17198d9e43cf7b
      obj:
        cmd: echo hello world && exit 1
        end: '2018-08-29T14:24:25.324076'
        meta: {}
        returncode: 1
        start:
        status: failed
        tid: b58e7bf6bc1cfff2437491baea17198d9e43cf7b
    - id: c53235a60e74239d95ab4927ba954399892c9782
      obj:
        cmd: echo hello world && exit 0
        end:
        meta: {}
        returncode:
        start:
        status: new
        tid: c53235a60e74239d95ab4927ba954399892c9782


How do we go about getting rid of that old task which we know will
always fail?

- A remove-failed command line tool that finds and removes all failed tasks?
- Always, or optionally, remove failed tasks from the workflow when starting?

Note that there may be cases where a task fails due to external circumstances, 
so retrying the task later is probably a common need.







