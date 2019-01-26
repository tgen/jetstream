import json
import logging
import os
import traceback
import jetstream

log = logging.getLogger(__name__)


class InvalidPipeline(Exception):
    pass


class ValidationFailed(Exception):
    pass


class Pipeline(object):
    """Represents a pipeline installed on the system. This will load the
    manifest and validate the data. To add a new validation, create a method
    with the name "validate_<field>" where field is the manifest value to
    validate. It should take one argument the value to validate."""
    sentinel = object()

    def __init__(self, path):
        self.path = os.path.abspath(path)
        self._manifest = Pipeline.sentinel
        self._constants = Pipeline.sentinel

    def __repr__(self):
        return f'<Pipeline({self.manifest.get("name")}) at {self.path}>'

    @property
    def manifest_path(self):
        filename = jetstream.settings['pipelines']['manifest_filename'].get()
        return os.path.join(self.path, filename)

    @property
    def constants_path(self):
        constants = jetstream.settings['pipelines']['constants_filename'].get()
        return os.path.join(self.path, constants)

    @property
    def constants(self):
        if self._constants is Pipeline.sentinel:
            self._constants = jetstream.load_file()
        return self._constants

    @property
    def manifest(self):
        """Provides lazy loading of Pipeline details"""
        if self._manifest is Pipeline.sentinel:
            self.reload()
        return self._manifest

    def load_manifest(self):
        manifest = self.manifest_path

        if not os.path.exists(manifest):
            raise InvalidPipeline(f'Manifest not found: {manifest}')

        with open(manifest, 'r') as fp:
            try:
                return json.load(fp)
            except json.decoder.JSONDecodeError:
                tb = traceback.format_exc()
                msg = f'Failed to load {self.path} error below:\n{tb}'
                raise InvalidPipeline(msg) from None

    def reload(self):
        if not os.path.isdir(self.path):
            raise InvalidPipeline(f'"{self.path}" is not a directory!')

        self._manifest = self.load_manifest()

        try:
            name = manifest['name']
        except KeyError:
            raise InvalidPipeline(f'manifest missing "name"') from None

        if not name.isidentifier():
            raise InvalidPipeline(f'"{name}" is not a valid identifier')


def ls():
    home = jetstream.settings['pipelines']['home']
    manifest = jetstream.settings['pipelines']['manifest_filename']

    pipelines = {}
    for d in os.listdir(home):
        path = os.path.join(home, d)
        manifest = os.path.join(path, manifest)

        if os.path.isdir(path) and os.path.exists(manifest):
            try:
                pipeline = Pipeline(path)
            except InvalidPipeline:
                log.exception('Failed to load pipeline!')

            pipelines[pipeline.name] = pipeline

    return pipelines


def get(name):
    pipelines = jetstream.pipelines.ls()

    try:
        pipeline = pipelines[name]
    except KeyError:
        manifest = {p.name: p.manifest for p in pipelines.values()}
        if manifest:
            msg = f'"{name}" not found.\n\nPipelines available:\n' \
                  f'{jetstream.utils.yaml_dumps(manifest)}'
        else:
            msg = 'No pipelines found!'

        raise ValueError(msg) from None

    if CONSTANTS.exists():
        return pipeline, CONSTANTS
    else:
        return pipelines, None


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
