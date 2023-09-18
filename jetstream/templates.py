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
import jinja2
from jinja2 import (
    Environment,
    StrictUndefined,
    Undefined,
    FileSystemLoader
)

log = logging.getLogger(__name__)

if hasattr(jinja2, "pass_context"):
    pass_context = jinja2.pass_context
else:
    pass_context = jinja2.contextfunction


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
            self.sources.append(f'Project: {rep}')
            self.stack.append(project.index)

        if pipeline is not None:
            ctx = pipeline.get_context()
            rep = textwrap.shorten(str(ctx), 76)
            self.sources.append(f'Pipeline: {rep}')
            self.stack.append(ctx)

        if command_args is not None:
            rep = textwrap.shorten(str(command_args), 76)
            self.sources.append(f'Command: {rep}')
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


@pass_context
def raise_helper(ctx, msg):
    """Allow "raise('msg')" to be used in templates"""
    raise TemplateException(f'{ctx.name}: {msg}')


@pass_context
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


def md5(path):
    """Allow "{{ path|md5 }}" to be used in templates. A good
    use case is to track the md5sum of a script or other file that may
    change over time. Causing the render to update on file change"""
    hash_md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def assignbin(value, bins=[0, float('inf')], labels=None):
    """Allow "{{ value|assignbin }}" to be used in templates. The return value
    of this filter is the 0-based bin the value falls in. Default bin is 0 to
    infinity. Edges floor to lower bin. Also accepts a list of labels such
    that: assignbin(5,[0,2,4,6],['low','med','high']) }} returns 'high'.
    Returns -1 if the value is out of bounds."""
    start = bins[:-1]
    end = bins[1:]
    for i, (lower_bound, upper_bound) in enumerate(zip(start, end)):
        if lower_bound <= value <= upper_bound:
            if labels is not None:
                return labels[i]
            return i
    return -1


def fromjson(value):
    """Allow "{{ value|fromjson }}" to be used in templates"""
    return json.loads(value)


def env(value):
    return os.environ[value]
    

def getenv(value, default=None):
    return os.environ.get(value, default)


def setenv(key, value):
    os.environ[key] = value
    return '' 


def environment(*searchpath, strict=True, trim_blocks=True, lstrip_blocks=True):
    """Starts a Jinja2 Environment with a FileSystemLoader on the given search
    path. This adds several features to the standard template processor."""
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
    env.globals['getenv'] = getenv
    env.globals['setenv'] = setenv
    env.filters['fromjson'] = fromjson
    env.filters['basename'] = basename
    env.filters['dirname'] = dirname
    env.filters['urlparse'] = urlparse
    env.filters['sha256'] = sha256
    env.filters['md5'] = md5
    env.filters['assignbin'] = assignbin
    return env


def load_template(path, *searchpath, **kwargs):
    """Helper function to quickly load a template from a file. Remaining args and kwargs
    are passed to jetstream.templates.environment() """
    template_dir = os.path.dirname(path)
    template_name = os.path.basename(path)
    env = environment(template_dir, *searchpath, **kwargs)
    return env.get_template(template_name)


def from_string(data, *searchpath, **kwargs):
    """Helper function to quickly load a template from a string. Remaining args and kwargs
    are passed to jetstream.templates.environment() """
    env = environment(*searchpath, **kwargs)
    return env.from_string(data)


def render_template(template, project=None, pipeline=None, command_args=None):
    """Render a template and return as a string"""
    log.info('Rendering template...')

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


def load_workflow(render):
    """Given a rendered template string, loads the tasks and returns a workflow"""
    log.debug(f'Parsing tasks from render:\n{render}')
    tasks = jetstream.utils.parse_yaml(render)

    if not tasks:
        raise ValueError('No tasks found!')

    log.info('Loading tasks...')
    if isinstance(tasks, Mapping):
        # If template is a mapping, it is expected to have a tasks section
        props = {k: v for k, v in tasks.items() if k != 'tasks'}
        tasks = [jetstream.Task(**t) for t in tasks['tasks']]
        wf = jetstream.Workflow(tasks=tasks, props=props)
    else:
        tasks = [jetstream.Task(**t) for t in tasks]
        wf = jetstream.Workflow(tasks=tasks)

    wf.reload_graph()
    return wf

