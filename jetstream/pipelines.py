import logging
import os
from distutils.version import LooseVersion
import jetstream

log = logging.getLogger(__name__)
MANIFEST_FILENAME = 'pipeline.yaml'


class InvalidPipeline(Exception):
    pass


class Pipeline(object):
    """Represents a pipeline installed on the system."""
    def __init__(self, path):
        self.path = os.path.abspath(os.path.expanduser(path))

        if not os.path.isdir(self.path):
            raise InvalidPipeline(f'"{self.path}" is not a directory')

        self.manifest_path = os.path.join(self.path, MANIFEST_FILENAME)
        self.manifest = jetstream.utils.load_file(self.manifest_path)

        try:
            self.name = self.manifest['name']
            self.version = str(self.manifest.get('version', 0))
        except KeyError as e:
            err = f'{self.path}: manifest missing key "{e}"'
            raise InvalidPipeline(err) from None

        if not self.name.isidentifier():
            err = f'{self.path}: "{self.name}" is not a valid identifier'
            raise InvalidPipeline(err)

    def __repr__(self):
        return f'<Pipeline: {self.name}-{self.version}>'


def find_pipelines(home=None):
    home = home or jetstream.settings['pipelines']['home'].get(str)
    for filename in os.listdir(home):
        path = os.path.join(home, filename)
        manifest_path = os.path.join(path, MANIFEST_FILENAME)
        if os.path.isdir(path) and os.path.exists(manifest_path):
            try:
                yield Pipeline(path)
            except Exception:
                log.exception(f'Failed to load: {path}')


def lookup(name, version=None, home=None):
    if version is not None:
        version = str(version)
        for p in find_pipelines(home):
            if p.name == name and p.version == version:
                return p
    else:
        matches = []
        for p in find_pipelines(home):
            if p.name == name:
                matches.append(p)

        if matches:
            s = sorted(matches, key=lambda p: LooseVersion(p.version))
            return s[-1]

    msg = f'Pipeline "{name}-{version}" not found!'
    raise FileNotFoundError(msg)


# TODO New feature - load workflows directly from a python module:
# # Tools for loading workflows from a python module
# def _find_modules(instance):
#     def is_mod(attr):
#         if attr.startswith('_'):
#             return False
#
#         try:
#             fn = getattr(instance, attr)
#         except AttributeError:
#             return False
#
#         if not callable(fn):
#             return False
#
#         return getattr(fn, 'is_module', False)
#
#     return filter(is_mod, dir(instance))
#
#
# def _pipeline_module(fn):
#     """Decorator applied to Pipeline methods to make them automatically
#     fire when the pipeline is instantiated"""
#     fn.is_module = True
#     return fn
#
#
# def _load_mod(filepath):
#     dirn = os.path.dirname(filepath)
#     base = os.path.basename(filepath)
#     name, ext = os.path.splitext(base)
#     old_path = sys.path.copy()
#     old_modules = sys.modules.copy()
#
#     try:
#         sys.path.insert(1, dirn)
#         spec = importlib.util.spec_from_file_location(name, filepath)
#         mod = importlib.util.module_from_spec(spec)
#         spec.loader.exec_module(mod)
#     finally:
#         sys.path = old_path
#         sys.modules = old_modules
#
#     return mod
#
#
# def _find_tasks(mod):
#     for k, v in mod.__dict__.items():
#         if isinstance(v, jetstream.Task):
#             yield v
