"""Initiate a Jinja environment with template loaders that search
locations set by arguments or environment variables. """
import os
import json
import hashlib
from types import MethodType
from jinja2 import (Environment, FileSystemLoader, meta,
                    StrictUndefined, Undefined, evalcontextfilter)
from jetstream import log



def raise_helper(msg):
    """Allow "raise('msg')" to be used in templates"""
    raise Exception(msg)


def basename(path):
    return os.path.basename(path)


def dirname(path):
    return os.path.dirname(path)


def sha256(value):
    h = hashlib.sha256(value.encode())
    return h.hexdigest()


@evalcontextfilter
def fromjson(eval_ctx, value):
    return json.loads(value)


def get_children(env, template):
    source, filename, loader = env.loader.get_source(env, template)
    parsed = env.parse(source)

    children = meta.find_referenced_templates(parsed)

    for c in children:
        yield c

        for d in get_children(env, c):
            yield d


def get_source(env, template):
    res = {template: env.loader.get_source(env, template)[:2]}

    for c in get_children(env, template):
        res[c] = env.loader.get_source(env, c)[:2]

    return res


def get_template_with_source(env, template, *args, **kwargs):
    t = env.get_template(template, *args, **kwargs)
    t.source = get_source(env, template)
    return t


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

    env = Environment(
        trim_blocks=trim_blocks,
        lstrip_blocks=lstrip_blocks,
        loader=FileSystemLoader(search_path),
        undefined=undefined_handler
    )

    env.get_template_with_source = MethodType(get_template_with_source, env)
    env.globals['raise'] = raise_helper
    env.filters['fromjson'] = fromjson
    env.filters['basename'] = basename
    env.filters['dirname'] = dirname
    env.filters['sha256'] = sha256

    log.debug('Template loader search path: {}'.format(env.loader.searchpath))

    return env


def load_template(path, template_search_path=None):
    # Configure the template environment and load
    template_name = os.path.basename(path)
    template_dir = os.path.dirname(path)
    search_path = [template_dir, ]

    if template_search_path:
        search_path.extend(template_search_path)

    env = environment(search_path=search_path)
    template = env.get_template_with_source(template_name)

    return template
