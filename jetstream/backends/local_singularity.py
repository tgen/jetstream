import os
import re
import glob
import shlex
import logging
import asyncio
import subprocess
import tempfile
import jetstream
from asyncio import Lock, BoundedSemaphore, create_subprocess_shell, CancelledError
from pathlib import Path

log = logging.getLogger('jetstream.backends')

class LocalSingularityBackend(jetstream.backends.BaseBackend):
    def __init__(self, cpus=None, blocking_io_penalty=None, pullfolder=None):
        """The LocalSingularityBackend executes tasks as processes on the local machine.

        This contains a semaphore that limits tasks by the number of cpus
        that they require. It requires that self.runner be set to get the
        event loop, so it's not instantiated until preflight.

        :param cpus: If this is None, the number of available CPUs will be
            guessed. This cannot be changed after starting the backend.
        :param blocking_io_penalty: Delay (in seconds) when a BlockingIOError
            prevents a new process from spawning.
        :param max_concurrency: Max concurrency limit
        """
        super(LocalSingularityBackend, self).__init__()
        self.cpus = cpus \
                    or jetstream.settings['backends']['local']['cpus'] \
                    or jetstream.utils.guess_local_cpus()
        self.bip = blocking_io_penalty \
                   or jetstream.settings['backends']['local']['blocking_io_penalty'].get(int)
        self._cpu_sem = BoundedSemaphore( self.cpus )
        memory_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        self.memory_gb = int( memory_bytes/(1024.**3) )
        self._mem_sem = BoundedSemaphore( self.memory_gb )
        self._resources_lock = Lock()
        self._singularity_run_sem = BoundedSemaphore( self.cpus ) # To ensure pulls have exclusive use of singularity
        self._singularity_pull_lock = Lock()
        self._singularity_pull_cache = {}
        log.info(f'LocalSingularityBackend initialized with {self.cpus} cpus and {self.memory_gb}G memory')
        
    async def spawn( self,
                     task,
                     allow_cpus_overbooking = True,
                     allow_memory_overbooking = True ):
        log.debug('Spawn: {}'.format(task))

        if 'cmd' not in task.directives:
            return task.complete()
    
        cmd = task.directives['cmd']
        input_filenames = task.directives.get( 'input', [] )
        output_filenames = task.directives.get( 'output', [] )
        cpus = task.directives.get('cpus', 0)
        
        memory_gb_required_str = task.directives.get('mem', "0G")
        memory_gb_required_value = int( memory_gb_required_str[:-1] )
        memory_gb_required_unit = memory_gb_required_str[-1]
        if memory_gb_required_unit == "M":
            memory_gb_required_value = 1
        elif memory_gb_required_unit != "G":
            raise RuntimeError('Task memory units must be M or G')

        docker_image = task.directives.get( 'docker_image', None )
        if docker_image == None:
            raise RuntimeError(f'docker_image argument missing for task: {task.name}')
        docker_image_split = docker_image.split()
        if len( docker_image_split ) > 1:
            singularity_image = " ".join(docker_image_split[:-1] + [ f"docker://{docker_image_split[-1]}" ] )
        else:
            singularity_image = f"docker://{docker_image}"
        
        cpus_reserved = 0
        memory_gb_reserved = 0
        open_fps = list()

        if cpus > self.cpus:
            if allow_cpus_overbooking:
                log.warning( f'Task cpus ({cpus}) greater than system cpus ({self.cpus}), {task.name}: Proceeding anyways ...' )
                cpus = self.cpus
            else:
                raise RuntimeError('Task cpus greater than system cpus')

        if memory_gb_required_value > self.memory_gb:
            if allow_memory_overbooking:
                log.warning( f'Task mem ({memory_gb_required_value}G) greater than system mem ({self.memory_gb}G), {task.name}: Proceeding anyways ...' )
                memory_gb_required_value = self.memory_gb
            else:
                raise RuntimeError('Task mem greater than system mem')
            
        log.debug( f'going to pull: {singularity_image}' )
        if singularity_image in self._singularity_pull_cache:
            pass
        else:
            async with self._singularity_pull_lock:
                for i in range( self.cpus ):
                    await self._singularity_run_sem.acquire()
                pull_command_run_string = f'singularity exec --cleanenv {singularity_image} true'
                log.debug( f'pulling: {pull_command_run_string}' )
                _p = await create_subprocess_shell( pull_command_run_string,
                                                    stdout=asyncio.subprocess.PIPE,
                                                    stderr=asyncio.subprocess.PIPE )
                stdout, stderr = await _p.communicate()
                log.debug( f'pulled, stdout: {stdout}' )
                log.debug( f'pulled, stderr: {stderr}' )
                self._singularity_pull_cache[ singularity_image ] = singularity_image
                for i in range( self.cpus ):
                    self._singularity_run_sem.release()
        log.debug( f'pulled: {singularity_image}' )
        
        try:
            async with self._resources_lock:
                log.debug('Reserving cpus: {}'.format(task))
                for i in range(cpus):
                    await self._cpu_sem.acquire()
                    cpus_reserved += 1
                
                log.debug('Reserving mem: {}'.format(task))
                for i in range(memory_gb_required_value):
                    await self._mem_sem.acquire()
                    memory_gb_reserved += 1
                
                log.debug('Resources reserved: {}'.format(task))
            
            stdin, stdout, stderr = self.get_fd_paths(task)

            if stdin:
                stdin_fp = open(stdin, 'r')
                open_fps.append(stdin_fp)
            else:
                stdin_fp = None

            if stdout:
                stdout_fp = open(stdout, 'w')
                open_fps.append(stdout_fp)
            else:
                stdout_fp = None

            if stderr:
                stderr_fp = open(stderr, 'w')
                open_fps.append(stderr_fp)
            else:
                stderr_fp = None
                    
            p = await self.subprocess_sh(
                cmd,
                task_name=task.name,
                input_filenames=input_filenames,
                output_filenames=output_filenames,
                singularity_image=singularity_image,
                stdin=stdin_fp,
                stdout=stdout_fp,
                stderr=stderr_fp
            )

            task.state.update(
                stdout_path=stdout,
                stderr_path=stderr,
                label=f'Slurm({p.pid})',
            )

            log.info(f'LocalSingularityBackend spawned({p.pid}): {task.name}')
            rc = await p.wait()

            if rc != 0:
                log.info(f'Failed: {task.name}')
                return task.fail(p.returncode)
            else:
                log.info(f'Complete: {task.name}')
                return task.complete(p.returncode)
        except CancelledError:
            task.state['err'] = 'Runner cancelled Backend.spawn()'
            return task.fail(-15)
        finally:
            for fp in open_fps:
                fp.close()

            for i in range(cpus_reserved):
                self._cpu_sem.release()

            for i in range(memory_gb_reserved):
                self._mem_sem.release()
                
            return task

    async def subprocess_sh( self, args, task_name, *, input_filenames=[], output_filenames=[],
                             singularity_image=None, stdin=None, stdout=None, stderr=None,
                             cwd=None, encoding=None, errors=None, env=None,
                             loop=None, executable='/bin/bash' ): 
        try:
            singularity_mounts = set()
            
            for input_filename_glob_pattern in input_filenames:
                input_filenames_glob = glob.glob( input_filename_glob_pattern )
                if len( input_filenames_glob ) == 0:
                    raise RuntimeError(f'input file(s) do not exist: {input_filename_glob_pattern}')
                for input_filename in input_filenames_glob:
                    input_filename = os.path.abspath( input_filename )
                    input_filename_head, input_filename_tail = os.path.split( input_filename )
                    singularity_mounts.add( input_filename_head )
            for output_filename in output_filenames:
                output_filename = os.path.abspath( output_filename )
                output_filename_head, output_filename_tail = os.path.split( output_filename )
                singularity_mounts.add( output_filename_head )
                Path( output_filename_head ).mkdir( parents=True, exist_ok=True )
            mount_strings = []
            for singularity_mount in singularity_mounts:
                mount_strings.append( "-B %s" % ( singularity_mount ) )
            singularity_mounts_string = " ".join( mount_strings )

            os.makedirs( "jetstream/cmd", mode = 0o777, exist_ok = True )
            run_script_filename = f"jetstream/cmd/{task_name}.bash"
            with open( run_script_filename, "w" ) as run_script:
                run_script.write( args )
            command_run_string = f"""singularity exec --cleanenv --nv {singularity_mounts_string} {self._singularity_pull_cache[ singularity_image ]} bash {run_script_filename}"""

            log.debug('command_run_string:\n------BEGIN------\n{}\n------END------'.format(command_run_string))
            
            while 1:
                try:
                    await self._singularity_run_sem.acquire()
                    log.debug('subprocess_run_sh: trying...')
                    p = await create_subprocess_shell(
                        command_run_string,
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
                    log.warning(f'System refusing new processes: {e}')
                    await asyncio.sleep(self.bip)
                finally:
                    self._singularity_run_sem.release()
                    
        except Exception as e:
            log.warning(f'Exception: {e}')
            p = await create_subprocess_shell(
                        "exit 1;",
                        stdin=stdin,
                        stdout=stdout,
                        stderr=stderr,
                        cwd=cwd,
                        encoding=encoding,
                        errors=errors,
                        env=env,
                        loop=loop,
                        executable=executable )
        return p
