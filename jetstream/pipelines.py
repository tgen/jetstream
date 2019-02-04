import logging
import os
import traceback
import yaml
import jetstream

log = logging.getLogger(__name__)
MANIFEST_FILENAME = 'pipeline.yaml'


class InvalidPipeline(Exception):
    pass


class Pipeline(object):
    """Represents a pipeline installed on the system."""
    def __init__(self, path):
        self.path = os.path.abspath(os.path.expanduser(path))
        self.manifest = self.load_manifest(self.path)

    def __repr__(self):
        if self.manifest is None:
            return f'<Pipeline at {self.path}>'
        else:
            return f'<Pipeline({self.manifest.get("name")}) at {self.path}>'

    @staticmethod
    def load_manifest(path):
        if not os.path.isdir(path):
            raise InvalidPipeline(f'"{path}" is not a directory!')

        manifest_path = os.path.join(path, MANIFEST_FILENAME)

        try:
            with open(manifest_path, 'r') as fp:
                manifest = jetstream.utils.yaml_loads(fp.read())
        except FileNotFoundError:
            err = f'Manifest file not found {manifest_path}'
            raise InvalidPipeline(err) from None
        except yaml.YAMLError:
            tb = traceback.format_exc()
            msg = f'Failed to load {path} error below:\n{tb}'
            raise InvalidPipeline(msg) from None

        try:
            name = manifest['name']
        except KeyError:
            raise InvalidPipeline(f'manifest missing "name"') from None

        if not name.isidentifier():
            raise InvalidPipeline(f'"{name}" is not a valid identifier')

        return manifest


def ls(home):
    pipelines = {}
    for d in os.listdir(home):
        path = os.path.join(home, d)
        manifest = os.path.join(path, MANIFEST_FILENAME)

        if os.path.isdir(path) and os.path.exists(manifest):
            try:
                pipeline = Pipeline(path)
            except InvalidPipeline:
                log.exception('Failed to load pipeline!')

            pipelines[pipeline.manifest['name']] = pipeline

    return pipelines


def lookup(name):
    home = jetstream.settings['pipelines']['home'].get()
    pipelines = jetstream.pipelines.ls(home)

    try:
        return pipelines[name]
    except KeyError:
        msg = f'Pipeline "{name}" not found!'
        raise ValueError(msg) from None


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
