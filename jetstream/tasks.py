from datetime import datetime
from hashlib import sha1
from jetstream import utils, log

valid_status = ('new', 'pending', 'complete', 'failed')


class InvalidTaskStatus(Exception):
    msg_prefix = "Invalid task status, options: {}: ".format(
        ', '.join(valid_status))

    def __init__(self, msg=''):
        super(InvalidTaskStatus, self).__init__(self.msg_prefix + msg)


class Task(object):
    def __init__(self, **kwargs):
        """Tasks are the fundamental units of a workflow.

        Tasks are composed of directives. Directives are key-value pairs given
        as kwargs when instantiating a new task. Directives are used by:

            Workflow - to build a graph of the task dependencies and ensure
            that tasks are executed in the correct order

            Backends - when executing tasks, directives control how the backend
            will reserve resources, record outputs, execute the command, etc.

        Tasks also have state attributes. Task.status is used by the workflow
        to coordinate tasks, while Task.state is meant to store additional data
        that should not be used to compare the content of two tasks.

        Task identity is computed by sha1 hash of the task directives, it can
        be accessed with Task.identity (note: this is not the same as the the
        object identity "id(Task)"). Task directives are immutable. Upon
        instantiation, they are stored as a tuple, json serialized, and hashed
        to generate an identity. The task directive tuple can be viewed with
        Task.identity. But, task directives should not be modified after
        instantiation, and attempting do so will have unpredictable effects on
        workflow execution. The ID is used to unambiguously compare tasks when
        building and combining Workflows.

        :param from_data: Rehydrate the task from existing data.
        :param kwargs: Generate a new task with the given directives.
        """
        self.workflow = None
        self._directives = tuple(kwargs.items())
        self._identity = hash_directives(self._directives)

        try:
            self.tid = kwargs['name']
        except KeyError:
            self.tid = self.identity

        self.clear_state()

    def __eq__(self, other):
        """Tasks can be compared to other Tasks, this also allows querying
        container objects with "in". Tasks are considered equal if their name
        is equal. When tasks have no name declared, their directives are used
        instead. """
        if isinstance(other, Task):
            return other.tid == self.tid
        else:
            return False

    def __hash__(self):
        return self.tid.__hash__()

    def __repr__(self):
        l = self.state.get("label")

        if l:
            return f'<Task({self.status}:{l}): {self.tid}>'
        else:
            return  f'<Task({self.status}): {self.tid}>'

    @property
    def status(self):
        return self.state['status']

    @status.setter
    def status(self, value):
        """Must be in Task.valid_status."""
        if value not in valid_status:
            raise InvalidTaskStatus(value)

        self.state['status'] = value

    @property
    def identity(self):
        """Hash of the directives that describes the task contents"""
        return self._identity

    def clear_state(self):
        log.debug(f'Clear state: {self}')
        directives = self.directives()

        # Retry is state-ish but provided in a directive. Each time the task is
        # instantiated, the remaining attempts are loaded into state again.
        # This means that retries will not persist across run instances.
        self.state = {
            'status': 'new',
            'remaining_attempts': directives.get('retry', 0)
        }

    def directives(self):
        """Access the directives as a dictionary, this is kinda slow because
        of the need to create a dictionary from the directives."""
        return dict(self._directives)

    def directives_raw(self):
        """Directives are stored as a tuple so that they're immutable"""
        return self._directives

    def reset(self, descendants=True, clear_state=True):
        """Reset the state of this task"""
        log.debug(f'Reset: {self}')
        if clear_state:
            self.clear_state()

        if self.workflow and descendants:
            for dep in self.descendants():
                dep.reset(descendants=False)

    def pending(self):
        """Indicate that this task has been passed to the runner"""
        log.debug(f'Pending: {self}')
        self.status = 'pending'
        self.state['start_time'] = datetime.now().isoformat()

    def fail(self, returncode=None, descendants=True):
        """Indicate that this task has failed"""
        log.debug(f'Failed: {self}')
        atts = self.state.get('remaining_attempts', 0)

        if atts > 0:
            self.reset(descendants=False)
            self.state['remaining_attempts'] = atts - 1
            return

        donedt = datetime.now()

        try:
            starts = self.state['start_time']
            startdt = datetime.strptime(starts,'%Y-%m-%dT%H:%M:%S.%f')
            elapsed = str(donedt - startdt)
        except (KeyError, ValueError):
            elapsed = None

        self.status = 'failed'
        self.state['returncode'] = returncode
        self.state['done_time'] = donedt.isoformat()
        self.state['elapsed_time'] = elapsed

        if self.workflow and descendants:
            for dep in self.descendants():
                dep.state.update(dependency_failed=self.tid)
                dep.fail(descendants=False)

    def complete(self, returncode=None):
        """Indicate that this task is complete"""
        log.debug(f'Complete: {self}')
        self.status = 'complete'
        self.state['done_at'] = datetime.now().isoformat()

        if returncode is not None:
            self.state.update(returncode=returncode)

    # Helpers and shortcut methods
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

    def serialize(self):
        return serialize(self)

    @staticmethod
    def deserialize(data):
        return deserialize(data)

    def copy(self):
        return copy(self)

    def flatten(self):
        return flatten(self)


def copy(task):
    return deserialize(serialize(task))


def deserialize(data):
    """Restores a task object that was previously exported as a dictionary with
    jetstream.tasks.serialize()"""
    task = Task(**data['directives'])
    task.state.update(data['state'])
    return task


def flatten(task):
    """Returns a flat dict representation of a task that may be useful for
    graph libraries that don't support complex node attributes."""
    res = {'tid': task.tid}

    for k, v in task.directives().items():
        res['directive_' + k] = v

    for k, v in task.state.items():
        res['state_' + k] = v

    return res


def hash_directives(directives):
    directives = utils.json_dumps(directives, sort_keys=True)
    return sha1(directives.encode('utf-8')).digest().hex()


def serialize(task):
    """Returns a task as a dictionary that can be easily dumped to YAML/json"""
    return {
        'tid': task.tid,
        'state': task.state,
        'directives': task.directives(),
    }
