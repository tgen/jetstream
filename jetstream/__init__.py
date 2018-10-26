import sys
import os
import re
import ulid
import traceback
from datetime import datetime
from jetstream import profile, logs
from jetstream.profile import load_profile, dumps_profile, profile_path

log = logs.logging.getLogger('jetstream')
settings = None


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


    def run_id(value=None):
        if value is not None:
            value = str(value)

            if value.isidentifier():
                return value
            else:
                raise ValueError('{} is not a valid identifier'.format(value))

        return run_id_template.format(ulid.new().str)


    def save_workflow(workflow, path):
        """Save a workflow to the path

        This helper function will try to choose the correct file format based
        on the extension of the path, but defaults to pickle for unrecognized
        extensions.

        :param workflow: Workflow instance
        :param path: where to save
        :return: None
        """
        format = workflow_extensions.get(os.path.splitext(path)[1], 'pickle')
        start = datetime.now()
        workflow_savers[format](workflow, path)
        elapsed = datetime.now() - start
        log.info('Workflow saved (after {}): {}'.format(elapsed, path))


    def load_workflow(path):
        """Load a workflow from a file.

        This helper function will try to choose the correct file format based
        on the extension of the path, but defaults to pickle for unrecognized
        extensions."""
        format = workflow_extensions.get(os.path.splitext(path)[1], 'pickle')
        return workflow_loaders[format](path)


    from jetstream import workflows, templates, projects, runner, backends
    from jetstream.runner import Runner
    from jetstream.projects import Project
    from jetstream.templates import (environment, render_template, load_template,
                                     render_templates, load_templates)
    from jetstream.workflows import (build_workflow, random_workflow, Workflow,
                                     Task)

    workflow_extensions = {
        '': 'pickle',
        'pickle': 'pickle',
        'yaml': 'yaml',
        'yml': 'yaml',
        'json': 'json',
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

