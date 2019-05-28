"""Initiate a Jinja2 environment with template loaders that search
locations set by arguments or environment variables. """
import json
import hashlib
import logging
import os
import urllib.parse
import textwrap
import jetstream
from jinja2 import (
    Environment,
    StrictUndefined,
    Undefined,
    contextfunction,
    FileSystemLoader
)

log = logging.getLogger(__name__)


@contextfunction
def raise_helper(ctx, msg):
    """Allow "raise('msg')" to be used in templates"""
    raise Exception(f'{ctx.name}: {msg}')


@contextfunction
def log_helper(ctx, msg):
    """Allow "raise('msg')" to be used in templates"""
    log.critical(f'{ctx.name}: {msg}')
    return ''

@contextfunction
def this_helper(ctx):
    return ctx


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
    env.globals['this'] = this_helper
    env.filters['fromjson'] = fromjson
    env.filters['basename'] = basename
    env.filters['dirname'] = dirname
    env.filters['urlparse'] = urlparse
    env.filters['sha256'] = sha256
    return env


class TemplateContext(object):
    """Enforces the correct priority order for loading config data sources when
    rendering templates:

    low priority  --> 1) Project: jetstream/config.yaml
                      2) Pipeline: pipeline.yaml: template_config_data section
    high priority --> 3) Command Args: -c/--config and -C/--config-file options

    """
    def __init__(self, *, project=None, pipeline=None, command_args=None):
       self.stack = (project.get_config, pipeline.config, command_args)

    def __repr__(self):
        sources = ' -> '.join((s for s in self.stack if s))
        return '<TemplateContext '

    def flatten(self):
        """Returns a single config object created by flattening the sources
        in """
        log.debug('Flattening template context stack...')
        stack = []

        if self.project:
            log.debug(f'Including project: {self.project.name}')
            stack.append(self.project.get_config())

        if self.pipeline:
            log.debug(f'Including pipeline: {self.pipeline.name}')
            stack.append(self.pipeline.get_config())

        if self.command_args:
            cmd_args = textwrap.shorten(str(self.command_args), 24)
            log.debug(f'Including cmd args: {cmd_args}')
            stack.append(self.command_args)

        return jetstream.utils.config_stack(*stack)


def render_template_file(path, env=None, **kwargs):
    """Using Jinja2 environment, loads a template from a filepath, then renders
    with the given context. See render_template for more info"""
    template_dir = os.path.dirname(path)
    template_basename = os.path.basename(path)

    if env is None:
        env = environment(searchpath=[template_dir, ])

    template = env.get_template(template_basename)
    return render_template(template, **kwargs)


def render_template_string(s, env=None, **kwargs):
    """Using Jinja2 environment, loads a template from a string, then renders
    with the given context. See render_template for more info"""
    if env is None:
        env = environment(searchpath=[os.getcwd(),])

    template = env.from_string(s)
    return render_template(template, **kwargs)


def render_template(template, project=None, pipeline=None, command_args=None):
    """Render a template with the given context data sources"""
    log.info('Rendering template...')
    context = TemplateContext(
        project=project,
        pipeline=pipeline,
        command_args=command_args
    )
    log.debug(f'Template render context: {context}')
    ctx = context.flatten()
    render = template.render(**ctx)
    return render
