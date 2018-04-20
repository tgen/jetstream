import os
import pkg_resources

__version__ = pkg_resources.get_distribution("jetstream").version

# This prevents numpy from starting a million threads when imported. The
# graph library, networkx, uses scipy/numpy. TODO switch to another graph lib
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

from jetstream import settings, utils, legacy
from jetstream.workflows import Workflow
from jetstream.project import Project


def load_project(path=None):
    return Project(path)


def load_workflow(path):
    wf_data = utils.yaml_load(path)
    return Workflow.from_node_link_data(wf_data)
