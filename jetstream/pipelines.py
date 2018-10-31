import os
import sys
import importlib
import jetstream
from jetstream import log


class Pipeline:
    def __init__(self):
        self.project = jetstream.Project()  # Variables specific to the current project
        self.constants = {}  # Variables constant all projects
        self.env = jetstream.environment()
        self.workflow = jetstream.Workflow()
        self.modules = []

    def build(self):
        log.info('Loading pipeline...')
        for mod in find_modules(self):
            getattr(self, mod)()

        for m in self.modules:
            for mod in find_modules(m):
                getattr(m, mod)(self)

        project_workflow = self.project.workflow()

        if project_workflow:
            project_workflow.compose(self.workflow)
            self.workflow = project_workflow
        else:
            self.workflow.save_path = self.project.workflow_file

        log.info('Pipeline loaded: {}'.format(self.workflow))

    def context(self):
        """The context is all variables for this run made available to render
        the directives for each task."""
        ctx = {
            'project': self.project,
            'config': self.project.config,
            'constants': self.constants
        }
        return ctx

    @staticmethod
    def _new_task(instance, context_locals=None, **directives):
        """Locals will be added to the context when rendering the task
        directives if locals is not None. This can be useful if you're
        looping through some variable adding tasks to the workflow, and
        the command needs to include some variable from the loop."""
        ctx = instance.context()

        if context_locals is not None:
            ctx.update(locals=context_locals)

        for k, v in directives.items():
            t = instance.env.from_string(v)
            r = t.render(ctx)
            directives[k] = r

        instance.workflow.new_task(**directives)

    def new_task(self, context_locals=None, **directives):
        self._new_task(self, context_locals=context_locals, **directives)


def find_modules(instance):
    def is_mod(attr):
        if attr.startswith('_'):
            return False

        try:
            fn = getattr(instance, attr)
        except AttributeError:
            return False

        if not callable(fn):
            return False

        return getattr(fn, 'is_module', False)

    return filter(is_mod, dir(instance))


def pipeline_module(fn):
    """Decorator applied to Pipeline methods to make them automatically
    fire when the pipeline is instantiated"""
    fn.is_module = True
    return fn


def load_mod(filepath):
    dirn = os.path.dirname(filepath)
    base = os.path.basename(filepath)
    name, ext = os.path.splitext(base)
    old_path = sys.path.copy()
    old_modules = sys.modules.copy()

    try:
        sys.path.insert(1, dirn)
        spec = importlib.util.spec_from_file_location(name, filepath)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.path = old_path
        sys.modules = old_modules

    return mod


def find_tasks(mod):
    for k, v in mod.__dict__.items():
        if isinstance(v, jetstream.Task):
            yield v


def load_pipeline(path):
    module = load_mod(path)
