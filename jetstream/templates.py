"""Initiate a Jinja environment with template loaders that search
locations set by arguments or environment variables. """
import os
import json
import logging
import jetstream
from types import MethodType
from jinja2 import (Environment, FileSystemLoader, meta,
                    StrictUndefined, Undefined, evalcontextfilter)

log = logging.getLogger(__name__)

SITE_TEMPLATES_ENVVAR = 'JETSTREAM_SITE_TEMPLATES'
STRICT_ENVVAR = 'JETSTREAM_STRICT'

SITE_TEMPLATES = bool(os.environ.get(SITE_TEMPLATES_ENVVAR, 1))
STRICT = bool(os.environ.get(STRICT_ENVVAR, 1))


def envbool(value):
    """Convert value to a boolean environment variable: str "0" or "1" """
    return str(int(bool(value)))


def raise_helper(msg):
    """Allow "raise('msg')" to be used in templates"""
    raise Exception(msg)


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


def get_template_with_source(self, template, *args, **kwargs):
    t = self.get_template(template, *args, **kwargs)
    t.source = get_source(self, template)
    return t


def load_environment(template_dirs=None, include_site_templates=SITE_TEMPLATES,
        strict=STRICT):
    """Start a Jinja2 Environment with the given template directories.

    Templates are loaded by a Jinja2 FilesystemLoader that includes built-in
    templates by default. The search path is also extended to include any
    directories in template_dirs."""
    search_path = list()

    if template_dirs:
        search_path.extend(template_dirs)

    if include_site_templates:
        search_path.append(jetstream.site_template_path)

    if strict:
        undefined_handler = StrictUndefined
    else:
        undefined_handler = Undefined

    os.environ[SITE_TEMPLATES_ENVVAR] = envbool(include_site_templates)
    os.environ[STRICT_ENVVAR] = envbool(strict)

    env = Environment(
        trim_blocks=True,
        lstrip_blocks=True,
        loader=FileSystemLoader(search_path),
        undefined=undefined_handler
    )

    env.get_template_with_source = MethodType(get_template_with_source, env)
    env.globals['raise'] = raise_helper
    env.filters['fromjson'] = fromjson

    log.debug('Template loader: {}'.format(env.loader.searchpath))
    return env
