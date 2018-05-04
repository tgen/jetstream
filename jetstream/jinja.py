"""Initiate a Jinja environment with template loaders that search
locations set by arguments or environment variables. """
import os
import json
from jinja2 import (Environment, PackageLoader, FileSystemLoader, ChoiceLoader,
    StrictUndefined, Undefined, evalcontextfilter)

project_loader = FileSystemLoader(os.path.join(os.getcwd(), 'templates'))
package_loader = PackageLoader('jetstream', 'templates')

PROJECT_TEMPLATES_ENVVAR = 'JETSTREAM_PROJECT_TEMPLATES'
STRICT_ENVVAR = 'JETSTREAM_STRICT'

PROJECT_TEMPLATES = bool(os.environ.get(PROJECT_TEMPLATES_ENVVAR, 1))
STRICT = bool(os.environ.get(STRICT_ENVVAR, 1))


def envbool(value):
    return str(int(bool(value)))


def raise_helper(msg):
    """Allows "raise('msg')" to be used in templates"""
    raise Exception(msg)


@evalcontextfilter
def fromjson(eval_ctx, value):
    return json.loads(value)


def template_env(include_project_templates=PROJECT_TEMPLATES, strict=STRICT):
    """Start a Jinja2 Environment with the given template directories.

    Templates are loaded by a Jinja2 ChoiceLoader that includes
    [<project>/templates, <package>/templates]. Project templates can
     be ignored by setting include_project_templates=False."""
    loaders = [package_loader]

    if include_project_templates:
        loaders.insert(0, project_loader)

    if strict:
        undefined_handler = StrictUndefined
    else:
        undefined_handler = Undefined

    os.environ[PROJECT_TEMPLATES_ENVVAR] = envbool(include_project_templates)
    os.environ[STRICT_ENVVAR] = envbool(strict)

    env = Environment(
        trim_blocks=True,
        lstrip_blocks=True,
        loader=ChoiceLoader(loaders),
        undefined=undefined_handler
    )

    env.globals['raise'] = raise_helper
    env.filters['fromjson'] = fromjson
    return env
