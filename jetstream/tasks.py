import logging
import re
from datetime import datetime
from hashlib import sha1
import jetstream

log = logging.getLogger(__name__)

# Identity should always include cmd and exec but may have other
# fields added later (ie container ids)
IDENTITY = (
    'cmd',
    'exec'
)

VALID_STATES = (
    'new',
    'pending',
    'complete',
    'failed',
    'skipped'
)


class TaskDirectiveProcessor:
    """Returns a callable that will preprocess known task directives. This
    reduces the guesswork that has to happen when using directives for workflow
    graphs or runner backends. Task directives are still open-ended, anything
    not listed in KNOWN will just get added to task.directives without any
    processing.

    TODO allow task backends to define known directives and check for collisons
    """
    KNOWN = {
        'cmd': 'none_or_str',
        'exec': 'none_or_str',
        'before': 'coerce_list',
        'after': 'coerce_list',
        'input': 'coerce_list',
        'output': 'coerce_list',
        'before-re': 'coerce_list',
        'after-re': 'coerce_list',
        'input-re': 'coerce_list',
        'reset': 'coerce_list',
    }

    def __call__(self, directives):
        for k, v in self.KNOWN.items():
            fn = self._get(v)
            directives[k] = fn(k, directives.get(k))

        return directives

    def _get(self, name):
        fn = getattr(self, name)

        if not callable(fn):
            err = f'DirectiveProcessor.{name} is not callable'
            raise ValueError(err)

        return fn

    def coerce_tuple(self, directive, value):
        return jetstream.utils.coerce_tuple(value)

    def coerce_list(self, directive, value):
        return jetstream.utils.coerce_list(value)

    def none_or_str(self,directive, value):
        if value is None:
            return None
        elif isinstance(value, str):
            return value
        else:
            err = f'{directive} directives must be str or None'
            raise ValueError(err)


class Task:
    processor = TaskDirectiveProcessor()
    name_regex = re.compile("^[A-Za-z0-9_-]*$")

    def __init__(self, name=None, state=None, **directives):
        self.directives = self.processor(directives)
        self.identity = self._get_identity(directives)
        self.state = state
        self.name = name or self.identity

        try:
            name_valid = self.name_regex.match(self.name)
        except Exception:
            raise ValueError(f'Could not validate task name: {self.name}')

        if not name_valid:
            err = f'Invalid task name "{self.name}". ' \
                  'Task names can only include alphanumeric, hypen, and ' \
                  'underscore characters'
            raise ValueError(err)

        if self.state is None:
            self.clear_state()

    def __eq__(self, other):
        """Tasks can be compared to other Tasks, this also allows querying
        container objects with "in". Tasks are considered equal if their names
        are equal. Note, a task can be equal by name, but different by identity
        """
        try:
            return other.name == self.name
        except AttributeError:
            return False

    def __hash__(self):
        """Allows this object to be used in hashed collections"""
        return self.name.__hash__()

    def __repr__(self):
        return f'<Task({self.status}): {self.name}>'

    def _get_identity(self, directives):
        components = []
        for d in IDENTITY:
            val = directives.get(d, '')

            try:
                val = str(val)
            except Exception:
                err = f'{d}: directive failed to string casting: {val}'
                raise ValueError(err)

            components.append(val)

        identity = ''.join(components)
        return sha1(identity.encode('utf8')).digest().hex()

    def _set_start_time(self):
        self.state['start_time'] = datetime.now().isoformat()

    def _set_done_time(self):
        done_dt = datetime.now()

        try:
            start = self.state.get('start_time', '')
            start_dt = datetime.strptime(start,'%Y-%m-%dT%H:%M:%S.%f')
            elapsed = str(done_dt - start_dt)
        except (KeyError, ValueError, TypeError):
            elapsed = None

        self.state['done_time'] = done_dt.isoformat()
        self.state['elapsed_time'] = elapsed

    @property
    def status(self):
        return self.state['status']

    @status.setter
    def status(self, value):
        """Setter for task status. Status is stored in the state dictionary
        and must be in VALID_STATES"""
        if value not in VALID_STATES:
            raise ValueError(f'Invalid task status: {value}')

        self.state['status'] = value

    def clear_state(self):
        """Reset the task state"""
        self.state = {
            'status': VALID_STATES[0],
            'remaining_attempts': self.directives.get('retry', 0)
        }

    def reset(self, clear_state=True):
        """Reset the state of this task"""
        log.debug(f'Reset: {self}')
        self.status = 'new'

        if clear_state:
            self.clear_state()

    def pending(self):
        """Set task status to pending
        Pending status indicates that the task has been passed to the runner
        and will be started soon. """
        log.debug(f'Pending: {self}')
        self.status = 'pending'
        self._set_start_time()

    def fail(self, returncode=None, force=False):
        """Set task status to failed

        :param returncode: Return code will be added to task.state
        :param descendants: Also fails any descendants of this task
        :param force: Ignore any remaining retry attempts
        :return:
        """
        log.debug(f'Failed: {self}')
        atts = self.state.get('remaining_attempts', 0)
        if atts and not force:
            self.reset()
            self.state['remaining_attempts'] = atts - 1
        else:
            self.status = 'failed'
            self.state['returncode'] = returncode
            self._set_done_time()

    def skip(self, **reasons):
        """Failed by the system due to dependencies or other reasons"""
        log.debug(f'Skipped: {self}')
        self.state.update(**reasons)
        self.status = 'skipped'

    def complete(self, returncode=None):
        """Indicate that this task is complete"""
        log.debug(f'Complete: {self}')
        self.status = 'complete'
        self._set_done_time()

        if returncode is not None:
            self.state['returncode'] = returncode

    def is_new(self):
        if self.status == 'new':
            return True
        else:
            return False

    def is_pending(self):
        if self.status == 'pending':
            return True
        else:
            return False

    def is_done(self):
        if self.status in ('complete', 'failed', 'skipped'):
            return True
        else:
            return False

    def is_complete(self):
        if self.status == 'complete':
            return True
        else:
            return False

    def is_failed(self):
        if self.status in ('failed', 'skipped'):
            return True
        else:
            return False

    def is_skipped(self):
        if self.status == 'skipped':
            return True
        else:
            return False

    def to_dict(self):
        return to_dict(self)

    def copy(self):
        return copy(self)


def copy(task):
    """Dirty hack to create a new task instance"""
    return from_dict(to_dict(task))

def random_task(name=None, input=None):
    if name is None:
        name = jetstream.guid(formatter='random_task_{id}')
    cmd = 'set -ue\n\nhostname\n\ndate\necho {} is running'
    output = name + '.txt'
    return Task(name=name, cmd=cmd, input=input, output=output)


def from_dict(data):
    """Restores a task object that was previously exported as a dictionary with
    jetstream.tasks.to_dict"""
    task = Task(**data)
    return task


def to_dict(task):
    """Returns a task as a dictionary that can be easily dumped to YAML/json."""
    return dict(name=task.name, state=task.state, **task.directives)
