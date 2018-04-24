import os
from jinja2 import (Environment, PackageLoader, FileSystemLoader, ChoiceLoader,
    StrictUndefined, Undefined)

project_loader = FileSystemLoader('templates')
package_loader = PackageLoader('jetstream_pipelines', 'templates')


def env(*dirs,
        include_project_templates=True,
        include_package_templates=True,
        strict=True):
    """Start a Jinja2 Environment with the given template directories.

    Templates are loaded via Jinja2 ChoiceLoader. This searches a list
    of loaders to find templates. In this case, any template directories
    given in arguments are preferred over the jetstream_piplines template
    directory."""

    paths = [os.path.realpath(d) for d in dirs]
    loaders = [FileSystemLoader(d) for d in paths]

    if include_project_templates:
        loaders.append(project_loader)

    if include_package_templates:
        loaders.append(package_loader)

    if strict:
        undefined_handler = StrictUndefined
    else:
        undefined_handler = Undefined

    return Environment(
        loader=ChoiceLoader(loaders),
        undefined=undefined_handler
    )
