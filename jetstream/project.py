import os
import ulid
import logging
from jetstream import utils

log = logging.getLogger(__name__)


class Project:
    def __init__(self, path):
        self.path = path
        self.project_data = os.path.join(self.path, 'project.json')
        self.run_data = os.path.join(self.path, '.jetstream')
        assert os.path.exists(self.project_data)
        assert os.path.exists(self.run_data)

    def new_run(self):
        run_id = ulid.new().str
        run_dir = os.path.join(self.run_data, run_id)
        log.critical('Initializing new run {}'.format(run_dir))
        os.mkdir(run_dir)

        run_data = {'created': utils.fingerprint()}

        run_data_path = os.path.join(run_dir, 'data')
        with open(run_data_path, 'w') as fp:
            utils.yaml.dump(run_data, fp, default_flow_style=False)

        return run_id


def find_latest_run(path=None):
    if path is None:
        path = os.getcwd()

    run_data_dir = os.path.join(path, '.jetstream')

    runs = []
    for d in os.listdir(run_data_dir):
        if _is_run(d):
            runs.append(d)
    return sorted(runs)[0]


def _is_run(path):
    """ Returns True if path is a valid run """
    if os.path.isdir(path) and _is_valid_run_id(path):
        return True
    else:
        return False


def _is_valid_run_id(id):
    """ Returns True if id is a valid run id """
    try:
        ulid.from_str(id)
        return True
    except (TypeError, ValueError):
        return False
