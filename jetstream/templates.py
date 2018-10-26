"""Initiate a Jinja environment with template loaders that search
locations set by arguments or environment variables. """
import os
import json
import hashlib
from jinja2 import (Environment, FileSystemLoader, StrictUndefined, Undefined,
                    evalcontextfilter)
from jetstream import log


@evalcontextfilter
def raise_helper(eval_ctx, msg):
    """Allow "raise('msg')" to be used in templates"""
    raise Exception(msg)


@evalcontextfilter
def basename(eval_ctx, path):
    """Allow "basename(<path>)" to be used in templates"""
    return os.path.basename(path)


@evalcontextfilter
def dirname(eval_ctx, path):
    """Allow "dirname(<path>)" to be used in templates"""
    return os.path.dirname(path)


@evalcontextfilter
def sha256(eval_ctx, value):
    """Allow "sha256(<value>)" to be used in templates"""
    h = hashlib.sha256(value.encode())
    return h.hexdigest()


@evalcontextfilter
def fromjson(eval_ctx, value):
    """Allow "fromjson(<value>)" to be used in templates"""
    return json.loads(value)


def environment(search_path=None, strict=True, trim_blocks=True,
                lstrip_blocks=True):
    """Start a Jinja2 Environment with the given template directories.

    Templates are loaded by a Jinja2 FilesystemLoader that includes built-in
    templates by default. The search path is also extended to include any
    directories in template_dirs."""
    if strict:
        undefined_handler = StrictUndefined
    else:
        undefined_handler = Undefined

    if search_path is None:
        search_path = [os.getcwd(),]

    env = Environment(
        trim_blocks=trim_blocks,
        lstrip_blocks=lstrip_blocks,
        loader=FileSystemLoader(search_path),
        undefined=undefined_handler
    )

    env.globals['raise'] = raise_helper
    env.filters['fromjson'] = fromjson
    env.filters['basename'] = basename
    env.filters['dirname'] = dirname
    env.filters['sha256'] = sha256

    log.info('Autoconfigured Jinja2 search path: {}'.format(env.loader.searchpath))
    return env


def load_template(path, search_path=None):
    """Load a template from a file path

    This will automatically configure a FileSystemLoader to search in the
    current directory and template directory.

    :param path: path to a template file
    :param search_path: a list of paths to add to search path
    :return: jinja2.Template
    """
    log.debug('Load template from: {}'.format(path))

    if os.path.isfile(path):
        log.debug('Template is a file')
    elif os.path.exists(path):
        log.debug('Template is non-file path that exists')

    template_name = os.path.basename(path)

    if search_path is None:
        search_path = [os.getcwd(), os.path.dirname(path)]

    env = environment(search_path=search_path)
    return env.get_template(template_name)


def render_template(path, data=None, search_path=None):
    """Load and render a template.

    :param data: Mapping of data used to render template
    :param args: Passed to load_template
    :param kwargs: Passed to load_template
    :return: Rendered template string
    """
    if data is None:
        data = dict()

    template = load_template(path, search_path=search_path)
    return template.render(**data)


def load_templates(*paths, search_path=None):
    """Load several templates, see load_template"""
    return [load_template(p, search_path=search_path) for p in paths]


def render_templates(*paths, data=None, search_path=None):
    """Load and render several templates, see render_template"""
    return [render_template(p, data=data, search_path=search_path)
            for p in paths]

