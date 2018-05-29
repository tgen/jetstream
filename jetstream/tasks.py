""" Tasks are intantiated by the runner when a workflow task is ready
to be executed. The runner keeps track of a queue of active task objects
while it's running a workflow. """
import os
import logging
from jetstream import utils

log = logging.getLogger(__name__)

SPEC = {
    'name': 'Jetstream task directives',
    'spec': {
        'task_properties': {
            'id': {
                'description': 'Unique id for the task',
                'required': True,
                'type': str
            },
            'type': {
                'description': 'Type of task to create (local, slurm, etc.)',
                'type': str
            },
            'cmd': {
                'description': 'Array of command arguments',
                'type': list
            },
            'stdin': {
                'description': 'Data to be sent to the command stdin',
                'type': str
            },
            'stdout': {
                'description': 'Write stdout to a file',
                'type': str
            },
            'stderr': {
                'description': 'Write stderr to a file',
                'type': str
            },
            'before': {
                'description': 'This task should execute before some other '
                               'task/s',
                'type': list,
            },
            'after': {
                'description': 'This task should execute after some other '
                               'task/s',
                'type': list,
            }
        },
    }
}


class BaseTask(object):
    def __init__(self, task_id, task_directives):
        self.task_id = task_id
        self.task_directives = task_directives
        self.stdout_path = None
        self.stderr_path = None
        self.stdin_data = None
        self.extras = dict()
        self.returncode = None
        self._proc = None

        if self.task_directives.get('stdout'):
            self.stdout_path = self.task_directives['stdout']
        else:
            default_filename = utils.cleanse_filename(self.task_id)
            self.stdout_path = os.path.join('logs', default_filename + '.log')

        self.stderr_path = self.task_directives.get('stderr', self.stdout_path)

        self.stdin_data = task_directives.get('stdin')

    @property
    def proc(self):
        if self._proc is None:
            raise AttributeError('This task has not been launched yet!')
        else:
            return self._proc

    @proc.setter
    def proc(self, value):
        self._proc = value

    def poll(self):
        raise NotImplementedError

    def kill(self):
        raise NotImplementedError

    def wait(self):
        raise NotImplementedError

    def launch(self):
        raise NotImplementedError


class SpecNameHere:
    def __init__(self, spec):
        self.spec = spec['spec']
        self.name = spec['name']
        self.task_properties = self.spec['task_properties']

    def required_props(self):
        for name, attrs in self.task_properties.items():
            if 'required' in attrs and attrs['required'] is True:
                yield name

    def has_unknown_props(self, task):
        errs = []
        for k in task.keys():
            if k not in self.task_properties:
                errs.append('Unknown task property "{}"'.format(k))
        return errs

    def missing_req_props(self, task):
        errs = []
        for prop in self.required_props():
            if prop not in task:
                errs.append('Missing required prop "{}"'.format(prop))
        return errs

    def wrong_prop_type(self, task):
        errs = []
        for k, v in task.items():
            if k not in self.task_properties:
                continue

            cls = self.task_properties[k]['type']
            if not isinstance(v, cls):
                errs.append('Property "{}" should be type "{}"'.format(k, cls))

        return errs

    def validate(self, task):
        errs = []
        errs += self.missing_req_props(task)
        errs += self.has_unknown_props(task)
        errs += self.wrong_prop_type(task)

        if errs:
            msg = '\n'.join(errs)
            log.critical('task validation failed for {}:\n{}'.format(task, msg))
            return False
        else:
            return True

    def coerce(self, task):
        # Stagein/out can be strings
        if 'stagein' in task:
            if isinstance(task['stagein'], str):
                task['stagein'] = {'src': task['stagein']}

        if 'stageout' in task:
            if isinstance(task['stageout'], str):
                task['stageout'] = {'src': task['stageout']}

        # Any task property that is supposed to be a list can
        # also be a string. The string will be coerced to a
        # len 1 list.
        for k, v in task.items():
            cls = self.task_properties[k]['_class']

            if cls is list and isinstance(v, str):

                log.warning('Coercing attribute "{}: {}" to list.'.format(
                    k, v, cls))
                v = [v, ]
            task[k] = cls(v)
