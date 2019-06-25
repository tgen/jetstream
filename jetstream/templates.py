"""Initiate a Jinja2 environment with template loaders that search
locations set by arguments or environment variables. """
import json
import hashlib
import logging
import os
import urllib.parse
import textwrap
from collections.abc import Mapping
import jetstream
from jinja2 import (
    Environment,
    StrictUndefined,
    Undefined,
    contextfunction,
    FileSystemLoader
)

log = logging.getLogger(__name__)


class TemplateContext:
    """Enforces the correct priority order for loading config data sources when
    rendering templates:

    low priority  --> 1) Project: jetstream/config.yaml
                      2) Pipeline: pipeline.yaml: template_config_data section
    high priority --> 3) Command Args: -c/--config and -C/--config-file options

    """
    def __init__(self, *, project=None, pipeline=None, command_args=None):
        self.sources = []
        self.stack = []

        if project is not None:
            rep = textwrap.shorten(str(project.index), 76)
            self.sources.append(f'Project index: {rep}')
            self.stack.append(project.index)

        if pipeline is not None:
            rep = textwrap.shorten(str(pipeline.manifest), 76)
            self.sources.append(f'Pipeline manifest: {rep}')
            self.stack.append(pipeline.manifest)

        if command_args is not None:
            rep = textwrap.shorten(str(command_args), 76)
            self.sources.append(f'Command args: {rep}')
            self.stack.append(command_args)

    def __str__(self):
        sources = [f'{i}: {text}' for i, text in enumerate(self.sources)]
        return '\n'.join(sources)

    def flatten(self):
        """Returns a single config object created by flattening the sources
        in """
        return jetstream.utils.config_stack(*self.stack)



class TemplateException(Exception):
    """Can be raised by the template itself using raise() function"""


@contextfunction
def raise_helper(ctx, msg):
    """Allow "raise('msg')" to be used in templates"""
    raise TemplateException(f'{ctx.name}: {msg}')


@contextfunction
def log_helper(ctx, msg, level='INFO'):
    """Allow "raise('msg')" to be used in templates"""
    level = logging._checkLevel(level)
    log.log(level, f'{ctx.name}: {msg}')
    return ''


def basename(path):
    """Allow "{{ path|basename }}" to be used in templates"""
    return os.path.basename(path)


def dirname(path):
    """Allow "{{ path|dirname }}" to be used in templates"""
    return os.path.dirname(path)


def urlparse(value):
    """Allow "{{ value|urlparse }}" to be used in templates. The return
    value of this filter is a named tuple that can be used to access parts
    of a qualified URL. eg: {{ path|urlparse|attr("netloc") }}"""
    return urllib.parse.urlparse(value)


def sha256(value):
    """Allow "{{ value|sha256 }}" to be used in templates"""
    h = hashlib.sha256(value.encode())
    return h.hexdigest()


def fromjson(value):
    """Allow "{{ value|fromjson }}" to be used in templates"""
    return json.loads(value)


def environment(strict=True, trim_blocks=True, lstrip_blocks=True,
                searchpath=None):
    """Starts a Jinja2 Environment with a FileSystemLoader on the given search
    path. This adds several features to the standard template processor.
    """
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
    env.globals['log'] = log_helper
    env.filters['fromjson'] = fromjson
    env.filters['basename'] = basename
    env.filters['dirname'] = dirname
    env.filters['urlparse'] = urlparse
    env.filters['sha256'] = sha256
    return env


def render_template(path=None, data=None, *, project=None, pipeline=None,
                   command_args=None, search_path=()):
    log.info('Rendering template...')

    if bool(path) == bool(data):
        raise TypeError('build template expected path or data argument')

    if path:
        template_dir = os.path.dirname(path)
        template_basename = os.path.basename(path)
        env = environment(searchpath=[template_dir,] + list(search_path))
        template = env.get_template(template_basename)
    else:
        env = environment(searchpath=list(search_path))
        template = env.from_string(data)

    context = TemplateContext(
        project=project,
        pipeline=pipeline,
        command_args=command_args
    )

    sources = textwrap.indent(str(context), " " * 4)
    if sources:
        log.info(f'Template rendering data sources include:\n{sources}')
    else:
        log.warning(f'No data for rendering template variables')

    context = context.flatten()
    return template.render(**context)



def build_template(*args, **kwargs):
    render = render_template(*args, **kwargs)
    tasks = jetstream.utils.yaml_loads(render)

    log.info('Loading tasks...')
    if isinstance(tasks, Mapping):
        # If template is a mapping, it is expected to have a tasks section
        props = {k: v for k, v in tasks.items() if k != 'tasks'}
        tasks = [jetstream.Task(**t) for t in tasks['tasks']]
        wf = jetstream.Workflow(tasks=tasks, props=props)
    else:
        tasks = [jetstream.Task(**t) for t in tasks]
        wf = jetstream.Workflow(tasks=tasks)

    return wf

