# Developers guide

New backends need to implement a coroutine method "spawn". They can optionally
implement a coroutine "coro" for managing additional background tasks. Here
is an example:

Note that you may want to handle the asyncio.CancelledError in these coroutines
otherwise an error will be reported when the runner exits early. 

```
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

```
