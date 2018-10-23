from datetime import datetime
from hashlib import sha1
from jetstream import utils, log


class Task(object):
    valid_status = ('new', 'pending', 'complete', 'failed')

    def __init__(self, *, from_data=None, **kwargs):
        """Tasks are the fundamental units of a workflow.

        Tasks are composed of directives. Directives are key-value pairs, given
        as kwargs when instantiating a new task. Directives are used by:

            Workflow - to build a graph of the task dependencies and ensure
            that tasks are executed in the correct order.

            Backends - when executing tasks, directives control how the backend
            will reserve resources, record outputs, execute the command, etc.

        Tasks also have state attributes. Task.status is used by the workflow
        to coordinate tasks, while Task.state is meant to store additional data
        during or after the task is run.

        Task identity is computed by sha1 hash of the task directives, it can
        be accessed with Task.tid (note: this is not the same as the the object
        identity "id(Task)"). Task directives are immutable. Upon instantiation
        they are stored as a tuple, then json serialized and hashed to generate
        an identity. The task directive tuple can be viewed with Task.identity.
        But, task directives should not be modified after instantiation, and
        attempting do so will have unpredictable effects on workflow execution.
        The ID is used to unambiguously compare tasks when building and
        executing workflows.

        :param from_data: Rehydrate the task from existing data.
        :param kwargs: Generate a new task with the given directives.
        """
        if from_data and kwargs:
            raise ValueError('from_data and kwargs are mutually exclusive')

        if from_data:
            self._tid = from_data['tid']
            self._state = utils.JsonDict(from_data['state'])
            self._directives = tuple(from_data['directives'].items())
            self.status = from_data['status']
        else:
            self._status = 'new'
            self._state = utils.JsonDict()
            self._directives = tuple(kwargs.items())
            identity = utils.json_dumps(self._directives, sort_keys=True)
            self._tid = sha1(identity.encode('utf-8')).digest().hex()

        self.workflow = None

    def __dict__(self):
        return {
            'tid': self.tid,
            'status': self.status,
            'directives': self.directives,
            'state': self.state.to_dict(),
        }

    def __eq__(self, other):
        if isinstance(other, Task):
            return other.tid == self.tid
        else:
            return other == self.tid

    def __hash__(self):
        return self.tid.__hash__()

    def __repr__(self):
        return '<Task {}>'.format(self.label)

    def serialize(self):
        return self.__dict__()

    def flat(self):
        """Returns a flat dict representation of this task that
        may be useful for graph libraries that don't support complex
        node attributes."""
        res = {'tid': self.tid, 'status': self.status}

        for k, v in self.directives.items():
            res['directive_' + k] = v

        for k, v in self.state.items():
            res['state_' + k] = v

        return res

    @property
    def tid(self):
        return self._tid

    @property
    def status(self):
        return self._status

    @property
    def label(self):
        name = self.directives.get('name')

        if name is not None:
            return name
        else:
            return self.tid

    @status.setter
    def status(self, value):
        """Must be in Task.valid_status."""
        if value not in Task.valid_status:
            raise InvalidTaskStatus(value)

        self._status = value

    @property
    def directives(self):
        """Access the directives as a dictionary"""
        return dict(self._directives)

    @property
    def identity(self):
        """Access the directives tuple"""
        return self._directives

    @property
    def state(self):
        """Access the state JsonDict"""
        return self._state

    def reset(self, quiet=False):
        """Reset the state of this task"""
        if not quiet:
            log.info('Reset: {}'.format(self))

        self.status = 'new'
        self._state = utils.JsonDict()

    def start(self, quiet=False):
        """Indicate that this task has been started"""
        if not quiet:
            log.info('Started: {}'.format(self))

        self.status = 'pending'
        self.state['start_at'] = datetime.now().isoformat()

    def fail(self, returncode=None, quiet=False):
        """Indicate that this task has failed"""
        if not quiet:
            log.info('Failed: {}'.format(self))

        self.status = 'failed'
        self.state['done_at'] = datetime.now().isoformat()

        if returncode is not None:
            self.state.update(returncode=returncode)

        if self.workflow:
            for dep in self.dependents():
                dep.state.update(dependency_failed=self.tid)
                dep.fail(quiet=True)

    def complete(self, returncode=None, quiet=False):
        """Indicate that this task is complete"""
        if not quiet:
            log.info('Complete: {}'.format(self))

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


class TaskException(Exception):
    """Base class for catching exceptions related to tasks"""
    msg_prefix = ''

    def __init__(self, msg=''):
        super(TaskException, self).__init__(self.msg_prefix + msg)


class InvalidTaskStatus(TaskException):
    msg_prefix = "Invalid task status, options: {}: ".format(
        ', '.join(Task.valid_status))
