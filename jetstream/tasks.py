import logging
from datetime import datetime
from hashlib import sha1
import jetstream

log = logging.getLogger(__name__)
valid_status = ('new', 'pending', 'complete', 'failed')


class InvalidTaskStatus(Exception):
    msg_prefix = f"Invalid task status, options: {', '.join(valid_status)}"

    def __init__(self, msg=''):
        super(InvalidTaskStatus, self).__init__(self.msg_prefix + msg)


class Task(object):
    def __init__(self, *, name=None, status='new', state=None, description=None,
                 **kwargs):
        """Tasks are the fundamental units of a workflow.

        Tasks are composed of directives. Directives are immutable key-value
        pairs given as kwargs when instantiating a new task. Directives are
        used by:

            Workflow - to build a graph of the task dependencies and ensure
            that tasks are executed in the correct order

            Runner/Backends - when executing tasks, directives control how the
            backend will reserve resources, save outputs, execute the
            command, etc.

        Tasks also have state attributes, state changes as the task runs.
        Task.state['status'] is used by the workflow to orchestrate workflow
        execution, it has a shortcut property Task.status. State does not
        influence the Task.identity.

        Task identity is computed by sha1 hash of the serialized task
        directives, it can be accessed with Task.identity (note: this is not
        the same as the the object identity "id(Task)"). Task directives are
        immutable. Upon instantiation, they are stored as a tuple, json
        serialized, and then hashed to generate the task identity. The task
        directive tuple can be viewed with Task.identity. But, task directives
        should not be modified after instantiation, and attempting do so will
        have unpredictable effects on workflow execution.

        :param kwargs: Task directives
        """
        self._directives = tuple(kwargs.items())
        self._identity = hash_directives(self._directives)
        self._status = None
        self.name = name or self.identity
        self.status = status
        self.state = state or dict()
        self.description = description or dict()
        self.workflow = None

    def __eq__(self, other):
        """Tasks can be compared to other Tasks, this also allows querying
        container objects with "in". Tasks are considered equal if their names
        are equal. """
        if isinstance(other, Task):
            return other.name == self.name
        else:
            return False

    def __hash__(self):
        """Allows this object to be used in hashed collections"""
        return self.name.__hash__()

    def __repr__(self):
        l = self.state.get("label")

        if l:
            return f'<Task({self.status}:{l}): {self.name}>'
        else:
            return f'<Task({self.status}): {self.name}>'

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
        return self._status

    @status.setter
    def status(self, value):
        """Must be in Task.valid_status."""
        if value not in valid_status:
            raise InvalidTaskStatus(value)

        self._status = value

    @property
    def identity(self):
        """Hash of the directives that describes the task contents"""
        return self._identity

    def clear_state(self):
        log.debug(f'Clear state: {self}')
        self.state = {}

    def directives(self):
        """Access the directives as a dictionary.

        This is intentionally left as a method, rather than a property, to
        emphasize the fact that it references a part of the task that should
        not be changed. """
        return dict(self._directives)

    def reset(self, descendants=True, clear_state=True):
        """Reset the state of this task"""
        log.debug(f'Reset: {self}')
        self.status = 'new'

        if clear_state:
            self.clear_state()

        if self.workflow and descendants:
            for dep in self.descendants():
                dep.reset(descendants=False)

    def pending(self):
        """Set task status to pending
        Pending status indicates that the task has been passed to the runner
        and will be started soon. """
        log.debug(f'Pending: {self}')
        self.status = 'pending'
        self._set_start_time()

    def fail(self, returncode=None, descendants=True, force=False):
        """Set task status to failed

        :param returncode: Return code will be added to task.state
        :param descendants: Also fails any descendants of this task
        :param force: Ignore any remaining retry attempts
        :return:
        """
        log.debug(f'Failed: {self}')

        if not force:
            if '_remaining_attempts' not in self.state:
                atts = self.directives().get('retry', 0)
                self.state['_remaining_attempts'] = atts

            if self.state['_remaining_attempts'] > 0:
                self.reset(descendants=False)
                self.state['_remaining_attempts'] -= 1
                return

        self.status = 'failed'
        self.state['returncode'] = returncode
        self._set_done_time()

        if self.workflow and descendants:
            for dep in self.descendants():
                dep.state.update(dependency_failed=self.name)
                dep.fail(descendants=False, force=True)

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
        if self.status in ('complete', 'failed'):
            return True
        else:
            return False

    def is_complete(self):
        if self.status == 'complete':
            return True
        else:
            return False

    def is_failed(self):
        if self.status == 'failed':
            return True
        else:
            return False

    def is_ready(self):
        return self.workflow.is_ready(self)

    def dependents(self):
        return self.workflow.dependents(self)

    def dependencies(self):
        return self.workflow.dependencies(self)

    def ancestors(self):
        return self.workflow.ancestors(self)

    def descendants(self):
        return self.workflow.descendants(self)

    def to_dict(self):
        return to_dict(self)

    @staticmethod
    def from_dict(data):
        return from_dict(data)

    def copy(self):
        return copy(self)


def copy(task):
    """Dirty hack to create a new task instance"""
    return from_dict(to_dict(task))


def from_dict(data):
    """Restores a task object that was previously exported as a dictionary with
    jetstream.tasks.to_dict"""
    task = Task(**data)
    return task


def hash_directives(directives):
    directives = jetstream.utils.json_dumps(directives, sort_keys=True)
    return sha1(directives.encode('utf-8')).digest().hex()


def random_task(name=None, input=None):
    if name is None:
        name = jetstream.guid(formatter='random_task_{id}')
    cmd = 'set -ue\n\nhostname\n\n\date\necho {} is running'
    output = name + '.txt'
    return Task(name=name, cmd=cmd, input=input, output=output)


def to_dict(task):
    """Returns a task as a dictionary that can be easily dumped to YAML/json"""
    res = task.directives()
    res.update({
        'name': task.name,
        'status': task.status,
        'state': task.state,
        'description': task.description,
    })

    return res
