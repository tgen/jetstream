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
    def __init__(self, *, project=None, pipeline=None, command_args=None, other_args=None):
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
        
        if other_args is not None:
            rep = textwrap.shorten(str(other_args), 76)
            self.sources.append(f'Arguments: {rep}')
            self.stack.append(other_args)

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


def cloudpath(path, dot_index=-1):
    parsed_path = urllib.parse.urlparse(path)
    if parsed_path.scheme:
        # This is a URL so pass it back unharmed
        return path
    
    # This is a directory path, so transform it into a cloud path
    split_path = path.rsplit(os.path.sep, dot_index * -1)
    split_path[dot_index] = './' + split_path[dot_index]
    return os.path.join(*split_path)


def cloudbasename(path, in_cloud):
    if in_cloud:
        return os.path.basename(path)
    return path


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
    env.filters['cloudpath'] = cloudpath
    env.filters['cloudbasename'] = cloudbasename
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


def render_template(template, project=None, pipeline=None, command_args=None, other_args=None):
    """Render a template and return as a string"""
    log.info('Rendering template...')

    context = TemplateContext(
        project=project,
        pipeline=pipeline,
        command_args=command_args,
        other_args=other_args
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

