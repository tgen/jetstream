import os
import pkg_resources

__version__ = pkg_resources.get_distribution("jetstream").version

project_index = 'jetstream'
project_config = 'config'
project_temp = 'temp'

# This prevents numpy from starting a million threads when imported. The
# graph library, networkx, uses scipy/numpy. TODO switch to another graph lib
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

from jetstream import settings, utils, legacy, workflows
from jetstream.workflows import Workflow
from jetstream.project import Project
from jetstream.jinja import template_env, package_loader


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


def load_project(path=None):
    return Project(path)


def load_workflow(path):
    wf_data = utils.yaml_load(path)
    return Workflow.from_node_link_data(wf_data)

