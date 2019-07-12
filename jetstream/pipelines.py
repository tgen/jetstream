import logging
import os
from distutils.version import LooseVersion
import jetstream

log = logging.getLogger(__name__)
MANIFEST_FILENAME = 'pipeline.yaml'


class InvalidPipeline(ValueError):
    pass


class Pipeline:
    """Represents a pipeline installed on the system.
    Pipelines are any directory with a manifest: pipeline.yaml. If the
    manifest does not contain required fields: "__pipeline__" with "name",
    "version" and "main", InvalidPipeline will be raised.

    Manifest file example:

    # This manifest is also available as config data
    # when rendering workflows. So, it can be a good
    # place to store extra variables that may be used
    # in the templates.

    __pipeline__:
      name: phoenix
      version: 1.0
      main: main.jst
      author: Ryan Richholt
      ...
    foo: bar
    constants:
      these: can be anything

    """
    def __init__(self, path):
        self.path = os.path.abspath(os.path.expanduser(path))

        if not os.path.isdir(self.path):
            raise InvalidPipeline(f'"{self.path}" is not a directory')

        self.manifest_path = os.path.join(self.path, MANIFEST_FILENAME)
        self.manifest = jetstream.utils.load_file(self.manifest_path)

        try:
            self.info = info = self.manifest['__pipeline__']
            self.name = info['name']
            self.main = info['main']
            self.version = str(info['version'])
        except KeyError as e:
            err = f'{self.path}: manifest missing key "{e}"'
            raise InvalidPipeline(err) from None

        if not self.name.isidentifier():
            err = f'{self.path}: "{self.name}" is not a valid identifier'
            raise InvalidPipeline(err)

    def __repr__(self):
        return f'<Pipeline: {self.name} ({self.version})>'


def is_pipeline(path):
    """Returns True if given path is likely a pipeline"""
    manifest_path = os.path.join(path, MANIFEST_FILENAME)
    if os.path.isdir(path) and os.path.exists(manifest_path):
        return True
    else:
        return False


def pipelines_iter(searchpath=None):
    """Yields all pipelines found in pipeline home (taken from settings file or
    defaulting to user home directory)."""
    settings_searchpath = jetstream.settings['pipelines']['searchpath'].get(str)
    searchpath = searchpath or settings_searchpath
    log.debug(f'Pipelines search path: {searchpath}')

    paths_to_search = searchpath.split(':')
    paths_to_search = [os.path.expanduser(p) for p in paths_to_search]

    for x in paths_to_search:
        log.debug(f'Searching for pipelines in: {x}')
        for filename in os.listdir(x):
            path = os.path.join(x, filename)
            if is_pipeline(path):
                try:
                    yield Pipeline(path)
                except Exception:
                    log.warning(f'Failed to load: {path}')


def list_pipelines(*args, **kwargs):
    """Returns all pipelines found in pipeline home as a list"""
    return list(pipelines_iter(*args, **kwargs))


def get_pipeline(name, version=None, *args, **kwargs):
    """Get a pipeline by name and version(optional)"""
    if version is not None:
        version = str(version)
        for p in pipelines_iter(*args, **kwargs):
            # TODO can we allow > < = syntax here?
            if p.name == name and p.version == version:
                return p
    else:
        matches = []
        for p in pipelines_iter(*args, **kwargs):
            if p.name == name:
                matches.append(p)

        if matches:
            s = sorted(matches, key=lambda p: LooseVersion(p.version))
            return s[-1]

    if version is not None:
        msg = f'Pipeline "{name} ({version})" not found!'
    else:
        msg = f'Pipeline "{name} (latest)" not found!'

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
