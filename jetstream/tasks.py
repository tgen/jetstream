from datetime import datetime
from hashlib import sha1
from jetstream import utils, log


class Task(object):
    valid_status = ('new', 'pending', 'complete', 'failed')
    def __init__(self, *, from_data=None, **kwargs):
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
        self._workflow = None

        if from_data and kwargs:
            raise ValueError('from_data and kwargs are mutually exclusive')

        if from_data:
            self._state = utils.JsonDict(from_data['state'])
            self._directives = tuple(from_data['directives'].items())
            self._directives_hash = self._hash_directives()
            self.status = from_data['status']
        else:
            self._state = utils.JsonDict()
            self._directives = tuple(kwargs.items())
            self._directives_hash = self._hash_directives()
            self.status = 'new'

        directives = self.directives()
        if 'name' in directives:
            self._tid = directives['name']
        else:
            self._tid = self._directives_hash

    def __eq__(self, other):
        """Tasks can be compared to other Tasks, this also allows querying
        container objects with "in"."""
        if isinstance(other, Task):
            return other.tid == self.tid
        else:
            return False

    def __hash__(self):
        return self.tid.__hash__()

    def __repr__(self):
        return f'<Task({self.status}): {self.tid}>'

    def _hash_directives(self):
        directives = utils.json_dumps(self._directives, sort_keys=True)
        return sha1(directives.encode('utf-8')).digest().hex()

    def copy(self):
        return Task.deserialize(self.serialize())

    def flat(self):
        """Returns a flat dict representation of this task that
        may be useful for graph libraries that don't support complex
        node attributes."""
        res = {'tid': self.tid, 'status': self.status}

        for k, v in self.directives().items():
            res['directive_' + k] = v

        for k, v in self.state.items():
            res['state_' + k] = v

        return res

    def serialize(self):
        return {
            'tid': self.tid,
            'status': self.status,
            'directives': self.directives(),
            'state': self.state.to_dict(),
        }

    @staticmethod
    def deserialize(from_data):
        return Task(from_data=from_data)

    @property
    def tid(self):
        return self._tid

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, value):
        """Must be in Task.valid_status."""
        if value not in Task.valid_status:
            raise InvalidTaskStatus(value)

        self._status = value

    @property
    def identity(self):
        """Hash of the directives"""
        return self._directives_hash

    @property
    def state(self):
        """Access the state JsonDict"""
        return self._state

    @property
    def workflow(self):
        return self._workflow

    @workflow.setter
    def workflow(self, value):
        if self._workflow is not None:
            raise AttributeError(
                f'{self} is already assigned to a workflow: {self._workflow}'
            )
        else:
            self._workflow = value

    def directives(self):
        """Access the directives as a dictionary, this is kinda slow because
        of the need to create a dictionary from the directives."""
        return dict(self._directives)

    def directives_raw(self):
        """Directives are stored as a tuple so that they're immutable"""
        return self._directives

    def reset(self, quiet=False, descendants=True):
        """Reset the state of this task"""
        if not quiet:
            log.info('Reset: {}'.format(self.tid))

        self.status = 'new'
        self._state = utils.JsonDict()

        if self.workflow and descendants:
            for dep in self.descendants():
                dep.reset(quiet=True)

    def pending(self, quiet=False):
        """Indicate that this task has been passed to the runner"""
        if not quiet:
            log.info('Pending: {}'.format(self.tid))

        self.status = 'pending'
        self.state['start_at'] = datetime.now().isoformat()

    def fail(self, returncode=None, quiet=False, descendants=True):
        """Indicate that this task has failed"""
        if not quiet:
            log.info('Failed: {}'.format(self.tid))

        self.status = 'failed'
        self.state['done_at'] = datetime.now().isoformat()

        if returncode is not None:
            self.state.update(returncode=returncode)

        if self.workflow and descendants:
            for dep in self.descendants():
                dep.state.update(dependency_failed=self.tid)
                dep.fail(quiet=True, descendants=False)

    def complete(self, returncode=None, quiet=False):
        """Indicate that this task is complete"""
        if not quiet:
            log.info('Complete: {}'.format(self.tid))

        self.status = 'complete'
        self.state['done_at'] = datetime.now().isoformat()

        if returncode is not None:
            self.state.update(returncode=returncode)

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


class TaskException(Exception):
    """Base class for catching exceptions related to tasks"""
    msg_prefix = ''

    def __init__(self, msg=''):
        super(TaskException, self).__init__(self.msg_prefix + msg)


class InvalidTaskStatus(TaskException):
    msg_prefix = "Invalid task status, options: {}: ".format(
        ', '.join(Task.valid_status))
