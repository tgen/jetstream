from copy import deepcopy
import os
import traceback
import textwrap
import jetstream
from jetstream import log, utils


kvarg_types = {
    "str": str,
    "int": int,
    "float": float,
    "bool": utils.to_bool,
    "json": utils.json_loads,
    "yaml": utils.yaml_loads,
    "file": jetstream.load_data_file,
    "file:csv": utils.load_table,
    "file:json": utils.load_json,
    "file:tab": utils.load_table,
    "file:table": utils.load_table,
    "file:tsv": utils.load_table,
    "file:yaml": utils.load_yaml,
}


class TemplateVariableArgumentError(Exception):
    pass


def _loader_helper(arg, argtype, key, value):
    """Loads template variables and gives informative error messages"""
    try:
        loader = kvarg_types[argtype]
    except KeyError:
        msg = f'Argument {arg}: Unknown arg type: {argtype}'
        raise TemplateVariableArgumentError(msg) from None

    try:
        obj = loader(value)
    except Exception:
        tb = traceback.format_exc()
        msg = f'Error details:\n\n' \
              f'{textwrap.indent(tb, "  ")}\n'\
              f'Variable loading failed for "{arg} {value}", traceback above' \

        raise TemplateVariableArgumentError(msg) from None

    return obj


def parse_kvargs(args, separator=':'):
    """Parse arbitrary key-value arguments: ``--<type>:<key> <value>``"""
    log.debug('Parsing kvargs: {}'.format(args))
    context = dict()
    iterator = iter(args)

    while 1:
        try:
            arg = next(iterator)
        except StopIteration:
            break

        if arg.startswith('--'):
            argtype, _, key = arg[2:].rpartition(separator)

            try:
                value = next(iterator)
            except StopIteration:
                msg = f'Argument {arg}: expected one argument'
                raise TemplateVariableArgumentError(msg) from None

            # "variables" is a special case for template variables keys that
            # updates the context directly instead of a specific attribute
            if key == 'variables':
                if argtype == '':
                    obj = _loader_helper(arg, 'file', key, value)
                else:
                    obj = _loader_helper(arg, argtype, key, value)

                try:
                    context.update(obj)
                except Exception:
                    msg = f'Error while updating render context with object ' \
                          f'loaded from "{arg} {value}". Loading "variables" ' \
                          f'should return a mapping object. '
                    raise TemplateVariableArgumentError(msg)

            elif key and argtype:
                context[key] = _loader_helper(arg, argtype, key, value)
            else:
                msg = f'Template variable argument parser encountered ' \
                      f'"{arg} {value}", which does not ' \
                      f'match the pattern for template variable arguments. ' \
                      f'Template variable args should be in the form: ' \
                      f'"--<type>:<key> <value>". If this argument was not ' \
                      f'meant to be a template variable, check for spelling ' \
                      f'errors or missing values, and review the command ' \
                      f'--help.'
                raise TemplateVariableArgumentError(msg)

    return context


def load_context(project=None, kvargs=None, kvarg_separator=None):
    if project:
        project = jetstream.Project(path=project)
        log.info(f'Working in {project}')
    else:
        try:
            project = jetstream.Project()
            log.info(f'Working in {project}')
        except jetstream.NotAProject:
            log.info('Not working in a project')
            project = None

    if project is None:
        context = parse_kvargs(kvargs, kvarg_separator)
        context['__project__'] = None
    else:
        project.config.update(parse_kvargs(kvargs, kvarg_separator))
        project.save_config()
        context = deepcopy(project.config)
        context['__project__'] = project

    return context
