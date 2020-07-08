from datetime import datetime
from abc import ABC, abstractmethod
import asyncio
from asyncio import create_subprocess_shell

class CloudStorageSession(ABC):
    def __init__(self):
        self._temp_container_name = 'jetstream-temp-{unique_datestamp}'.format(
            unique_datestamp=datetime.now().strftime('%Y%m%d%H%M%S%f')
        )
        self._container_created = False
    
    # @abstractmethod
    # def provider_login(self, *args):
    #     raise NotImplementedError
    
    @abstractmethod
    async def upload_blob(self, filepath, blobpath=None, container=None):
        raise NotImplementedError
    
    @abstractmethod
    def download_blob(self, filepath, blobpath=None, container=None):
        raise NotImplementedError
    
    async def subprocess_sh(
            self, args, *, stdin=None, stdout=None, stderr=None,
            cwd=None, encoding=None, errors=None, env=None,
            loop=None, executable="/bin/bash", blocking_io_penalty=10):
        """Asynchronous version of subprocess.run

        This will always use a shell to launch the subprocess, and it prefers
        /bin/bash (can be changed via arguments)"""
        # print('world')
        # log.debug(f'subprocess_sh:\n{args}')
        # log.info(f'subprocess_sh:\n{args}')
        # print('here')
        while True:
            try:
                p = await create_subprocess_shell(
                    args,
                    stdin=stdin,
                    stdout=stdout,
                    stderr=stderr,
                    cwd=cwd,
                    encoding=encoding,
                    errors=errors,
                    env=env,
                    loop=loop,
                    executable=executable
                )
                break
            except BlockingIOError as e:
                # log.warning(f'System refusing new processes: {e}')
                print(f'System refusing new processes: {e}')
                await asyncio.sleep(blocking_io_penalty)
            except Exception as e:
                print(f'Exception: {e}')