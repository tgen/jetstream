import logging
import os
from packaging.version import parse as parse_version
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

    Pipeline features:
        - Quickly reference workflow templates by a name and version
        - Stored config data is available when rendering the template
        - Stored supporting scripts/binaries are added to PATH when running
        - Stored environment variables are added to env when running

        TODO FUTURE:
        -  Packaged and deployed to remote hosts for backends without
           a shared filesystem. Otherwise all the included files would
           not be accessible when the commands run.

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
      env:
        FOO: BAZ
      ...
    foo: bar
    constants:
      these: can be anything

    """
    def __init__(self, path=None, validate=True):
        if path is None:
            path = os.getcwd()

        self.path = os.path.abspath(os.path.expanduser(path))
        self.name = None
        self.version = None
        self.main = None
        self.manifest = None
        self.env = None

        if validate:
            self.validate()

    def __repr__(self):
        if self.name is None:
            return f'<Pipeline(Not loaded yet): {self.path}>'
        return f'<Pipeline {self.name} ({self.version}): {self.path}>'

    def get_context(self):
        ctx = self.manifest.copy()
        ctx['__pipeline__']['path'] = self.path
        return ctx

    def load_manifest(self):
        manifest_path = os.path.join(self.path, MANIFEST_FILENAME)
        return jetstream.utils.load_file(manifest_path)

    def load_template(self):
        path = os.path.join(self.path, self.main)
        return jetstream.templates.load_template(path)

    def validate(self):
        """Loads the manifest, sets instance attributes, and validates the contents"""
        manifest = self.load_manifest()

        try:
            pipeline_info = manifest['__pipeline__']
        except KeyError as e:
            err = 'Pipeline manifest missing "__pipeline__" section'
            raise InvalidPipeline(err) from e

        try:
            name = pipeline_info['name']
        except KeyError as e:
            err = 'Pipeline info section missing "name" field'
            raise InvalidPipeline(err) from e

        if not name.isidentifier():
            err = f'{self.path}: "{name}" is not a valid identifier'
            raise InvalidPipeline(err)

        try:
            main = pipeline_info['main']
        except KeyError:
            main = 'main.jst'

        main_path = os.path.join(self.path, main)
        if not os.path.exists(main_path):
            err = f'Pipeline main not found: {main_path}'
            raise InvalidPipeline(err)

        try:
            version = str(pipeline_info['version'])
        except KeyError:
            version = '0.0.0'

        self.name = name
        self.version = version
        self.main = main
        self.manifest = manifest

    def set_environment_variables(self):
        if self.manifest is None:
            e = "Pipeline.set_environment_variables() called before validate()"
            raise ValueError(e)

        os.environ['JS_PIPELINE_PATH'] = self.path
        os.environ['JS_PIPELINE_NAME'] = self.name
        os.environ['JS_PIPELINE_VERSION'] = self.version

        bin_path = self.manifest['__pipeline__'].get('bin')
        if bin_path:
            bin_path = os.path.join(self.path, bin_path)
            new_path = f'{bin_path}:{os.environ["PATH"]}'
            os.environ['PATH'] = new_path
            os.environ['SINGULARITYENV_APPEND_PATH'] = bin_path

        if self.env:
            for k, v in self.env.items():
                os.environ[k] = v



def find_pipelines(*dirs):
    """Generator yields all pipelines in given directories. Any pipeline found will also
    be searched for nested pipelines. """
    for dirname in dirs:
        dirname = os.path.abspath(os.path.expanduser(dirname))

        if not os.path.isdir(dirname):
            log.debug(f'Not a directory that exists, skipping {dirname}')
            continue

        log.debug(f'Searching for pipelines in: {dirname}')
        for filename in os.listdir(dirname):
            path = os.path.join(dirname, filename)

            # Here we check if the item is likely to be a pipeline, this filters
            # out standard files, and prevents trying to instantiate and validate
            # a pipeline object for every directory we encounter.
            if is_pipeline(path):
                try:
                    p = Pipeline(path)
                    log.debug(f'Found {p} at {path}')
                    yield p
                except Exception:
                    log.debug(f'Failed to load: {path}')

                yield from find_pipelines(path)


def get_pipeline(name, version=None, searchpath=None):
    """Get a pipeline by name and version(optional)"""
    if searchpath is None:
        searchpath = jetstream.settings['pipelines']['searchpath'].get()
        searchpath = searchpath.split(':')

    if version is None:
        # Find all but then sort by version and return the latest
        matches = []
        for p in find_pipelines(*searchpath):
            if p.name == name:
                matches.append(p)

        if matches:
            s = sorted(matches, key=lambda p: parse_version(p.version))
            return s[-1]
    else:
        # Find a match with name and version
        for p in find_pipelines(*searchpath):
            # TODO can we allow > < = syntax here?
            if p.name == name and parse_version(str(p.version)) == parse_version(str(version)):
                return p


    if version is None:
        msg = f'Pipeline "{name}" not found!'
    else:
        msg = f'Pipeline "{name}({version})" not found!'

    raise FileNotFoundError(msg)


def is_pipeline(path):
    """Returns True if given path is likely a pipeline"""
    manifest_path = os.path.join(path, MANIFEST_FILENAME)
    if os.path.isdir(path) and os.path.exists(manifest_path):
        return True
    else:
        return False


def list_pipelines(*dirs):
    """Returns all pipelines found as a list"""
    return list(find_pipelines(*dirs))

