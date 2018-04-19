"""Jetstream is a collection of tools for automating workflows at TGen."""
import os

__version__ = '0.1.0a1'
__author__ = 'Ryan Richholt'
# TODO finish these attributes

# This prevents numpy from starting a million threads when imported. The
# graph library, networkx, uses scipy/numpy. TODO switch to another graph lib
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

from jetstream import settings, utils
from jetstream.core.workflows.workflow import Workflow
from jetstream.core import project
from jetstream.core.project import Project


def load_project(path=None):
    return Project(path)


def load_workflow(path):
    wf_data = utils.yaml_load(path)
    return Workflow.from_node_link_data(wf_data)
