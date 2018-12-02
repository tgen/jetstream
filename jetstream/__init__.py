import sys
import os
import re
import ulid
import traceback
from jetstream import profile, logs
from jetstream.profile import load_profile, dumps_profile, profile_path

__name__ = 'jetstream'
log = logs.logging.getLogger(__name__)
settings = None
_current_project = None


def run_id():
    return run_id_template.format(ulid.new().str)


try:
    settings = load_profile()

    # Configure parallel library dependencies (Used by numpy)
    if 'OPENBLAS_NUM_THREADS' not in os.environ:
        _openblas_num_threads = str(settings['openblas_num_threads'])
        os.environ.update(OPENBLAS_NUM_THREADS=_openblas_num_threads)

    if 'MKL_NUM_THREADS' not in os.environ:
        _mkl_num_threads = str(settings['mkl_num_threads'])
        os.environ.update(MKL_NUM_THREADS=_mkl_num_threads)

    # Sortable unique id generator and precompiled regex pattern
    run_id_template = 'js{}'
    run_id_pattern = re.compile(r'^js[A-Z0-9]{26}$')

    from jetstream import utils
    from jetstream.utils import load_data_file, loadable_files
    from jetstream import workflows, templates, projects, runner, backends
    from jetstream.runner import Runner
    from jetstream.projects import Project, NotAProject
    from jetstream.pipelines import Pipeline, pipeline_module
    from jetstream.workflows import (Workflow, Task, load_workflow,
                                     save_workflow,
                                     random_workflow,)
    from jetstream.templates import environment, render_template

    workflow_extensions = {
        '': 'pickle',
        '.pickle': 'pickle',
        '.yaml': 'yaml',
        '.yml': 'yaml',
        '.json': 'json',
    }

    workflow_loaders = {
        'pickle': workflows.load_workflow_pickle,
        'yaml': workflows.load_workflow_yaml,
        'json': workflows.load_workflow_json
    }

    workflow_savers = {
        'pickle': workflows.save_workflow_pickle,
        'yaml': workflows.save_workflow_yaml,
        'json': workflows.save_workflow_json
    }


except Exception as e:
    msg = 'Current settings profile:\n{}\n'\
          'Error! Jetstream failed to load, there may be errors in the ' \
          'settings profile.\nSee traceback for details.'
    print(traceback.format_exc(), file=sys.stderr)
    print(msg.format(dumps_profile(settings)), file=sys.stderr)
    sys.exit(1)

