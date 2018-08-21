from jetstream import utils
from datetime import datetime
from hashlib import sha1


class Task(object):
    state_values = ('status', 'returncode', 'start', 'end')
    valid_status = ('new', 'pending', 'complete', 'failed')

    def __init__(self, data=None, **kwargs):
        """
        TODO: Note in docs, access to task data is restricted. Changes to task
        state should be handled through methods, and changing task directives
        might have unpredictable consequences to their workflow.

        :param data: Load task from pre-existing data
        :param kwargs: New task with 
        """
        self.workflow = None
        self._directives = utils.JsonDict()
        self._state = utils.JsonDict()
        
        existing_id = None
        
        if data:
            # Status is the most important attribute of a task. An invalid
            # status could break the workflow iterator, so we double check 
            # the status of tasks that are being reloaded from data.
            if 'status' in data and not data['status'] in Task.valid_status:
                raise InvalidTaskStatus(data['status'])
            
            for k, v in data.items():
                if k == 'id':
                    existing_id = v
                elif k in Task.state_values:
                    self._state[k] = v
                else:
                    self._directives[k] = v
        else:
            self._directives.update(kwargs)
            self.reset()

        self._tid = sha1(self.identity.encode()).digest().hex()

        if existing_id and existing_id != self._tid:
            raise ValueError(
                'Different id from rehydrated task:\n{}\n{}\n'
                'Task data: {}'.format(existing_id, self._tid, data
                                       ))

    def __repr__(self):
        name = self.get('name', self.tid)
        return '<Task({}):{}>'.format(name, self.status)

    def __hash__(self):
        return hash(self._tid)

    def __getitem__(self, item):
        return self._directives.__getitem__(item)
    
    def __setitem__(self, key, value):
        raise ValueError('Task data should not be modified.')
    
    def __eq__(self, other):
        if hasattr(other, 'tid'):
            return self.tid == other.tid
        else:
            return False

    def get(self, *args, **kwargs):
        return self._directives.get(*args, **kwargs)

    @property
    def identity(self):
        return utils.json_dumps(self._directives, sort_keys=True)
 
    @property
    def tid(self):
        # Task.tid is deliberately not "Task.id" in order to try and prevent
        # confusing the "id(Task)", the object identity, with "Task.tid", the 
        # content identity computed by sha1 hash of Task.identity 
        return self._tid

    def to_json(self):
        res = dict(tid=self.tid)
        res.update(self._directives)
        res.update(self._state)
        return utils.json_dumps(res, sort_keys=True)

    @property
    def state(self):
        return self._state
    
    @property
    def status(self):
        return self._state.get('status', 'new')

    @status.setter
    def status(self, value):
        if value not in Task.valid_status:
            raise InvalidTaskStatus(value)

        self._state['status'] = value

    # Methods to update task state
    def reset(self):
        self._state.update(
            status='new',
            returncode=None,
            start=None,
            end=None
        )

    def start(self):
        """Indicate that this task has been started"""
        self._state.update(
            status='pending',
            start=datetime.now().isoformat()
        )

    def complete(self, returncode=0):
        """Indicate that this task is complete"""
        self._state.update(
            status='complete',
            returncode=returncode,
            end=datetime.now().isoformat(),
        )
        
    def fail(self, returncode=1):
        self._state.update(
            status='failed',
            returncode=returncode,
            end=datetime.now().isoformat()
        )
        
        if self.workflow:
            for dep in self.workflow.dependents(self):
                dep.set_state(dependency_failed=str(self))
                dep.fail(returncode=123)

    def set_state(self, **kwargs):
        """Add information to task state meta data"""
        if 'meta' not in self._state:
            self._state['meta'] = utils.JsonDict()

        self._state['meta'].update(**kwargs)

    # Helper methods for checking status
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


class TaskException(Exception):
    """Base class for catching exceptions related to tasks"""
    pass


class InvalidTaskStatus(TaskException):
    def __init__(self, msg):
        msg = "Invalid task status: '{}'. Options: {}".format(
            msg, ', '.join(Task.valid_status))
        super(InvalidTaskStatus, self).__init__(msg)
