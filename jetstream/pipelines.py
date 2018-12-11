import os
import sys
import json
import importlib
import traceback
from pathlib import Path
import jetstream
from jetstream import log

MANIFEST_FILENAME = 'pipeline.json'
CONSTANTS_FILENAME = 'constants.yaml'

try:
    JETSTREAM_PIPELINES = Path(os.environ['JETSTREAM_PIPELINES']).resolve()
    CONSTANTS = JETSTREAM_PIPELINES.joinpath(CONSTANTS_FILENAME)
except KeyError:
    JETSTREAM_PIPELINES = None
    CONSTANTS = None

class InvalidPipeline(Exception):
    pass


class ValidationFailed(Exception):
    pass


class Pipeline(object):
    """Represents a pipeline installed on the system. This will load the
    manifest and validate the data. To add a new validation, create a method
    with the name "validate_<field>" where field is the manifest value to
    validate. It should take one argument the value to validate."""
    def __init__(self, path):
        self.path = Path(path).resolve()

        if not self.path.is_dir():
            raise InvalidPipeline(f'"{path}" is not a directory!')

        self.manifest_file = self.path.joinpath(MANIFEST_FILENAME)

        if self.manifest_file.is_file():
            with open(self.manifest_file, 'r') as fp:
                try:
                    self.manifest = json.load(fp)
                except json.decoder.JSONDecodeError:
                    msg = f'Failed to load {self.path} error below:'
                    log.exception(msg)
                    raise InvalidPipeline from None
        else:
            raise InvalidPipeline(f'Manifest not found: {MANIFEST_FILENAME}')

        self.name = self.manifest['name']

        if not self.name.isidentifier():
            raise ValidationFailed(f'"{self.name}" is not a valid identifier')

        self.main = self.path.joinpath(self.manifest['main'])

        if 'bin' in self.manifest:
            self.bin = self.path.joinpath(self.manifest['bin'])
        else:
            self.bin = None

    def __repr__(self):
        return f'<Pipeline({self.manifest.get("name")}) at {self.path}>'


def ls():
    if JETSTREAM_PIPELINES is None:
        raise ValueError('PIPELINES_DIR is not set!')

    pipelines = {}
    for d in os.listdir(str(JETSTREAM_PIPELINES)):
        path = os.path.join(JETSTREAM_PIPELINES, d)
        manifest = os.path.join(path, MANIFEST_FILENAME)

        if os.path.isdir(path) and os.path.exists(manifest):
            try:
                pipeline = Pipeline(path)
            except InvalidPipeline:
                pass

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


# Tools for loading workflows from a python module
def _find_modules(instance):
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


def _pipeline_module(fn):
    """Decorator applied to Pipeline methods to make them automatically
    fire when the pipeline is instantiated"""
    fn.is_module = True
    return fn


def _load_mod(filepath):
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


def _find_tasks(mod):
    for k, v in mod.__dict__.items():
        if isinstance(v, jetstream.Task):
            yield v
