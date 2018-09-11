from datetime import datetime
from hashlib import sha1
from copy import deepcopy
from jetstream import utils, log


class Task(object):
    state_values = ('status', 'returncode', 'start', 'end', 'meta')
    valid_status = ('new', 'pending', 'complete', 'failed')

    def __init__(self, data=None, **kwargs):
        """Tasks are the fundamental unit of a workflow

        Task.tid is deliberately not "Task.id" in order to prevent confusion.
        The Python object identity "id(Task)" is not the same as "Task.tid".
        Task.tid is the identity computed by sha1 hash the task directive
        content. The raw data used to generate the hash can be viewed with
        Task.identity, and the directives used to generate the data can be
        viewed with Task.directives. But, these are just copied of the data
        and modifying the actual data should be avoided.

        Changes to task state should be handled through methods, and changing
        task directives might have unpredictable consequences to their workflow.

        :param data: Load task from pre-existing data
        :param kwargs: New task with 
        """
        self.workflow = None
        self._directives = utils.JsonDict()
        self._state = utils.JsonDict()
        self._tid = ''
        existing_tid = None
        
        if data:
            log.debug('Rehydrating task from: {}'.format(data))

            if 'status' in data and not data['status'] in Task.valid_status:
                raise InvalidTaskStatus(data['status'])

            for k, v in data.items():
                if k == 'tid':
                    # Consume the existing tid to check it against the new tid
                    existing_tid = v
                elif k == 'meta':
                    self._state['meta'] = utils.JsonDict(v)
                elif k in Task.state_values:
                    log.debug('Set state key: "{}" value: "{}"'.format(k, v))
                    self._state[k] = v
                else:
                    log.debug('Set directive key: "{}" value: "{}"'.format(k, v))
                    self._directives[k] = v
        else:
            self._directives.update(kwargs)
            self._init_state()

        self._set_tid()

        if existing_tid and existing_tid != self._tid:
            msg = (
                'Conflicting tid from rehydrated task!\n'
                'tid from data: {}\n'
                'self tid: {}\n'
                'self identity: {}\n'
                'self directives: {}\n'
                'self state: {}\n'
            )
            raise ValueError(msg.format(
                existing_tid,
                self.tid,
                self.identity,
                self.directives,
                self.state
            ))

    def __repr__(self):
        return '<Task {}>'.format(self.tid[:8])

    def __hash__(self):
        return hash(self._tid)

    def __eq__(self, other):
        if hasattr(other, 'tid'):
            return self.tid == other.tid
        else:
            return False

    def _set_tid(self):
        self._tid = sha1(self.identity.encode()).digest().hex()

    def _init_state(self):
        self._state.update(
            status='new',
            returncode=None,
            start=None,
            end=None,
            meta=utils.JsonDict()
        )

    def reset(self):
        """Reset the state state of this task"""
        log.info('{} is reset'.format(self))
        self._init_state()

    def start(self):
        """Indicate that this task has been started"""
        log.info('{} has started'.format(self))
        self._state.update(
            status='pending',
            start=datetime.now().isoformat()
        )

    def complete(self, returncode=0):
        """Indicate that this task is complete"""
        log.info('{} is complete'.format(self))
        self._state.update(
            status='complete',
            returncode=int(returncode),
            end=datetime.now().isoformat(),
        )

    def fail(self, returncode=1):
        log.info('{} has failed'.format(self))
        self._state.update(
            status='failed',
            returncode=int(returncode),
            end=datetime.now().isoformat()
        )

        if self.workflow:
            for dep in self.workflow.dependents(self):
                dep.set_state(dependency_failed=self.tid)
                dep.fail(returncode=123)

    def done(self, returncode):
        rc = int(returncode)
        if rc != 0:
            self.fail(rc)
        else:
            self.complete(rc)

    def set_state(self, **kwargs):
        """Add information to task state meta data"""
        log.info('{} set state: {}'.format(self, kwargs))
        if 'meta' not in self._state:
            self._state['meta'] = utils.JsonDict()

        self._state['meta'].update(**kwargs)

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

    @property
    def directives(self):
        return deepcopy(self._directives)

    @property
    def state(self):
        return deepcopy(self._state)

    @property
    def identity(self):
        return utils.json_dumps(self.directives, sort_keys=True)
 
    @property
    def tid(self):
        return self._tid

    @property
    def status(self):
        return self._state.get('status', 'new')

    @status.setter
    def status(self, value):
        if value not in Task.valid_status:
            raise InvalidTaskStatus(value)

        self._state['status'] = value

    def serialize(self):
        res = dict(tid=self.tid)
        res.update(self.directives)

        state = self.state
        if 'meta' in state:
            state['meta'] = state['meta'].to_dict()
        res.update(state)

        return res


class TaskException(Exception):
    """Base class for catching exceptions related to tasks"""
    msg_prefix = ''

    def __init__(self, msg=''):
        super(TaskException, self).__init__(self.msg_prefix + msg)


class InvalidTaskStatus(TaskException):
    msg_prefix = "Invalid task status, options: {}: ".format(
        ', '.join(Task.valid_status))
