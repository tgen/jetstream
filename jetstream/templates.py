"""Initiate a Jinja environment with template loaders that search
locations set by arguments or environment variables. """
import json
import hashlib
import logging
import os
import traceback
import confuse
from datetime import datetime
import jetstream
from jinja2 import (
    Environment,
    StrictUndefined,
    Undefined,
    evalcontextfilter,
    FileSystemLoader
)

log = logging.getLogger(__name__)


@evalcontextfilter
def raise_helper(eval_ctx, msg):
    """Allow "raise('msg')" to be used in templates"""
    raise Exception(msg)


@evalcontextfilter
def basename(eval_ctx, path):
    """Allow "basename(<path>)" to be used in templates"""
    if isinstance(path, str):
        return os.path.basename(path)
    else:
        return [os.path.basename(p) for p in path]


@evalcontextfilter
def dirname(eval_ctx, path):
    """Allow "dirname(<path>)" to be used in templates"""
    # TODO Make this work on sequences of paths
    if isinstance(path, str):
        return os.path.dirname(path)
    else:
        return [os.path.dirname(p) for p in path]


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


def context(*, project=None, command_args=None, settings=True):
    """Enforces the correct priority order for loading data sources when
    rendering templates"""

    # Template rendering config stack
    #
    #
    # constants:
    #   command args "constants" key
    #   pipeline manifest "constants" section (if using pipelines)
    #   settings file "constants" section
    #
    # <globals>: <values>
    #   command args
    #   project config file (if working in a project)
    #
    #
    #Constants load order: Config files > Pipeline manifest (pipelines only) > command args???"""]
    # Load constants
    stack = []

    if settings:
        sc = jetstream.settings.flatten().get('constants')
        if sc:
            stack.append(sc)

    if command_args:
        cc = command_args.get('constants', {})
        if cc:
            stack.append(cc)

    a = {'constants': jetstream.utils.config_stack(*stack)}

    # Load globals
    stack = []

    if project:
        stack.append(project.config)

    if command_args:
        stack.append(command_args)

    b = jetstream.utils.config_stack(*stack)

    # Mash them all together
    return jetstream.utils.config_stack(a, b)


def render_template(path, context=None, env=None, render_only=False):
    """Load and render a template.

    :param context: Mapping object with data used to render template
    :param env: A preconfigured Jinja Environment
    :param render_only: Return the rendered template instead of a workflow
    :return: jetstream.Workflow
    """
    log.info('Rendering template...')
    started = datetime.now()

    if env is None:
        env = environment(searchpath=[os.getcwd(), os.path.dirname(path)])

    if not os.path.exists(path):
        raise FileNotFoundError(path)

    if os.path.isdir(path):
        raise FileNotFoundError(f'{path} is a directory!')

    with open(path, 'r') as fp:
        data = fp.read()

    log.debug(f'Template render context:\n{context}')
    template = env.from_string(data)
    render = template.render(**context)

    if render_only:
        return render

    log.debug(f'Rendered template:\n{render}')
    log.info('Parsing template...')

    try:
        parsed_tasks = jetstream.utils.yaml_loads(render)
    except jetstream.utils.yaml.YAMLError:
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

    # log.debug(f'Parsed template:\n{parsed_tasks}')
    log.info('Building workflow...')

    workflow = jetstream.Workflow()
    with workflow:
        for task in parsed_tasks:
            if not isinstance(task, dict):
                msg = f'Error with task:\n{task}' \
                      f'The template was rendered and parsed successfully, ' \
                      f'but encountered a task that was not a mapping type.'
                raise ValueError(msg)

            workflow.new_task(**task)

    elapsed = datetime.now() - started
    log.info(f'Workflow ready (after {elapsed}): {workflow}')
    return workflow
