"""Initiate a Jinja environment with template loaders that search
locations set by arguments or environment variables. """
import os
import json
import hashlib
import traceback
from datetime import datetime
from jinja2 import (
    Environment,
    StrictUndefined,
    Undefined,
    evalcontextfilter,
    FileSystemLoader
)
from jetstream import log, utils
from jetstream.workflows import Workflow


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
    # TODO Make this work on sequences of paths
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


def environment(strict=True, trim_blocks=True, lstrip_blocks=True,
                searchpath=None):
    """Start a Jinja2 Environment with the given template directories.

    Templates are loaded by a Jinja2 FilesystemLoader that includes built-in
    templates by default. The search path is also extended to include any
    directories in template_dirs."""
    if strict:
        undefined_handler = StrictUndefined
    else:
        undefined_handler = Undefined

    if searchpath is None:
        searchpath = [os.getcwd(),]

    env = Environment(
        trim_blocks=trim_blocks,
        lstrip_blocks=lstrip_blocks,
        undefined=undefined_handler,
        loader=FileSystemLoader(searchpath=searchpath),
        extensions=['jinja2.ext.do']
    )

    env.globals['raise'] = raise_helper
    env.filters['fromjson'] = fromjson
    env.filters['basename'] = basename
    env.filters['dirname'] = dirname
    env.filters['sha256'] = sha256
    return env


def render_template(path, variables=None, env=None):
    """Load and render a template.

    :param variables: Mapping of data used to render template
    :param env: A preconfigured Jinja Environment
    :return: Rendered template string
    """
    log.info('Building workflow...')
    started = datetime.now()

    if variables is None:
        variables = dict()

    if env is None:
        env = environment(searchpath=[os.getcwd(), os.path.dirname(path)])

    if not os.path.exists(path):
        raise FileNotFoundError(path)

    if os.path.isdir(path):
        raise FileNotFoundError(f'{path} is a directory!')

    with open(path, 'r') as fp:
        data = fp.read()

    template = env.from_string(data)
    render = template.render(**variables)

    log.debug(f'Rendered Template:\n{render}')

    try:
        parsed_tasks = utils.yaml_loads(render)
    except utils.yaml.YAMLError:
        log.critical(f'Error parsing rendered template:\n{render}')
        tb = traceback.format_exc()
        log.critical(f'Parser Traceback:\n{tb}\n')
        msg = (
            f'The template was rendered successfully, but failed to parse with'
            f'the YAML loader. See the parser traceback above for more details.'
        )
        raise ValueError(msg) from None

    if not isinstance(parsed_tasks, list):
        log.warning(f'Error loading rendered template:\n{render}')
        msg = (
            f'The template was rendered and parsed successfully, but the '
            f'resulting data structure was a {type(parsed_tasks)} . Templates '
            f'should always produce a list of tasks.'
        )
        raise ValueError(msg)

    workflow = Workflow()

    with workflow:
        for task in parsed_tasks:
            if not isinstance(task, dict):
                log.critical(f'Error with task:\n{task}')
                msg = (
                    f'The template was rendered and parsed successfully, but '
                    f'encountered a task that was not a mapping type.'
                )
                raise ValueError(msg)

            workflow.new_task(**task)

    elapsed = datetime.now() - started
    log.info('Workflow ready after {}'.format(elapsed))
    return workflow
