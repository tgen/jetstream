import os
import logging
from pkg_resources import get_distribution, resource_filename

__version__ = get_distribution('jetstream').version

VERBOSE = 5
logging.addLevelName(VERBOSE, 'VERBOSE')

def verbose(self, message, *args, **kws):
    if self.isEnabledFor(VERBOSE):
        self._log(VERBOSE, message, args, **kws)

logging.Logger.verbose = verbose

built_in_templates = resource_filename('jetstream', 'built_in_templates')
run_id_template = 'js{}'
project_index = 'jetstream'
project_config = 'config'
project_temp = 'temp'
project_logs = 'logs'
project_pid_file = os.path.join(project_index, 'pid')
project_manifest = os.path.join(project_index, 'manifest')
project_workflow = os.path.join(project_index, 'workflow')
project_history = os.path.join(project_index, 'history')

# This prevents numpy from starting a bunch of threads when imported. The
# graph library, networkx, uses scipy/numpy. TODO switch to another graph lib?
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

log = logging.getLogger(__name__)
log.setLevel(1)

from jetstream.exc import *
from jetstream import utils, legacy


data_loaders = {
    '.txt': utils.table_to_records,
    '.csv': utils.table_to_records,
    '.mer': utils.table_to_records,
    '.tsv': utils.table_to_records,
    '.json': utils.json_load,
    '.yaml': utils.yaml_load,
    '.yml': utils.yaml_load,
    '.config': legacy.config.load,
}


from jetstream.runner import AsyncRunner, SlurmBackend, LocalBackend
from jetstream.workflows import Workflow, Task
from jetstream.projects import Project
from jetstream import templates, projects, workflows

project_init = Project.init
load_template = templates.load_template
build_workflow = workflows.build_workflow
template_environment = templates.environment


def load_data_file(path):
    """Attempts to load a data file from path, raises :ValueError
    if an suitable loader function is not found in data_loaders"""
    for ext, fn in data_loaders.items():
        if path.endswith(ext):
            loader = fn
            break
    else:
        raise ValueError('No loader fn found for {}'.format(path))

    return loader(path)


def loadable_files(directory):
    """Generator yields all files we can load (see data_loaders) """
    for file in os.listdir(directory):
        path = os.path.join(directory, file)
        if os.path.isfile(path) \
                and path.endswith(tuple(data_loaders.keys())):
            yield path
