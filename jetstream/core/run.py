"""Run records get saved every time a workflow or plugin is launched in a
project. This module contains functions for examining those records."""
import os

import ulid
from jetstream import exc, utils
from jetstream.core import settings


class Run(object):
    def __init__(self, path):
        self.path = path

    @property
    def logs(self):
        for path in self._logs():
            yield {path: utils.yaml_load(path)}
        return

    def _logs(self):
        log_path = os.path.join(self.path, 'logs')
        for dirpath, _, filenames in os.walk(log_path):
            for f in filenames:
                yield os.path.join(log_path, dirpath, f)
        return

    @property
    def workflow(self):
        try:
            wf = utils.yaml_load(self._workflow())
        except FileNotFoundError:
            return dict()
        return wf

    def _workflow(self):
        return os.path.join(self.path, 'workflow')

    @property
    def created(self):
        try:
            c = utils.yaml_load(self._created())
        except FileNotFoundError:
            return dict()
        return c

    def _created(self):
        return os.path.join(self.path, 'created')

    @property
    def data(self):
        obj = dict()
        obj.update(self.created)
        obj.update(self.workflow)
        obj.update({'logs': [log for log in self.logs]})
        return obj


def load_run(path):
    """Helper function for loading a run, validates the path before
    loading."""
    run_id = os.path.basename(path)
    if is_valid_run_id(run_id):
        return Run(path)
    else:
        raise exc.NotARun


def new_run_id():
    return settings.profile['RUN_DATA_PREFIX'] + ulid.new().str


def is_run(path):
    """ Returns True if path is a valid run """
    if os.path.isdir(path) and is_valid_run_id(os.path.basename(path)):
        return True
    else:
        return False


def is_valid_run_id(run_id):
    """ Returns True if id is a valid run id """
    prefix = settings.profile['RUN_DATA_PREFIX']
    try:
        if run_id.startswith(prefix):
            ulid.from_str(run_id[len(prefix):])
            return True
        else:
            return False
    except (TypeError, ValueError):
        return False
