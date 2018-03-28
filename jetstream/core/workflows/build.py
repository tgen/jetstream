import yaml
import json
import logging
from jinja2 import Template, Environment, meta, StrictUndefined, Undefined
from jetstream.core.workflows import spec
from jetstream.core.workflows.workflow import Workflow


log = logging.getLogger(__name__)


class NodeValidationError(Exception):
    pass


class WorkflowBuilderError(Exception):
    pass


def validate_nodes(nodes, coerce=True):
    errors = []
    if coerce:
        log.critical('Automatic type coercion enabled')
        for node in nodes:
            spec.current.coerce(node)

    for node in nodes:
        valid = spec.current.validate(node)
        if not valid:
            errors.append(node['id'])

    if errors:
        msgs = ', '.join([str(e) for e in errors])
        raise NodeValidationError('Validation failed for:\n{}'.format(msgs))
    else:
        return True


def render_template(template, obj=None, strict=False):
    """Use Jinja2 to render a yaml template with a given object"""
    if obj is None:
        obj = dict()

    # Warn about missing vars before we even try to render
    env = Environment()
    ast = env.parse(template)
    variables = meta.find_undeclared_variables(ast)
    for v in variables:
        if not v in obj:
            log.warning('Undeclared template variable: {}'.format(v))

    # Allow strict rendering
    if strict:
        t = Template(template, undefined=StrictUndefined)
    else:
        t = Template(template, undefined=Undefined)

    res = t.render(**obj)
    return res


def parse_nodes(nodes):
    """Given a sequence of nodes (dictionaries with properties described in
     the workflow specification), returns a workflow with nodes and edges
     built. """
    wf = Workflow()

    for node in nodes:
        log.debug('Adding node to workflow: {}'.format(node))
        wf.add_node(**node)

    for node in nodes:
        # Before/After linking
        if 'before' in node:
            log.debug('Adding dependency to node: {}'.format(node))
            wf.add_dependency(node['id'], before=node['before'])

        if 'after' in node:
            wf.add_dependency(node['id'], after=node['after'])

        # In/Out linking
        # if 'in' in node:
        #     for dep in node['in']:
        #         found_match = False
        #         for n in wf.nodes(data=True):
        #             if n['out'] == dep:
        #                 wf.add_dependency(node['id'], after=(n,))
        #                 found_match = True
        #         if found_match is False:
        #             msg = 'Unmatched "in": {}'.format(node)
        #             raise WorkflowBuilderError(msg)

        # TODO This type of linking is much more complicated
        # need to revisit and think about how to address edge
        # cases, ex out directive with no ins, number of matches
        # per in/out etc.

    return wf


def build(template, data):
    # Render template with project_data using Jinja2
    log.critical('Rendering template...')
    raw = render_template(template, data)

    # Load template with yaml
    log.critical('Parsing template...')
    nodes = yaml.load(raw)

    # validate nodes
    log.critical('Validating nodes...\n{}'.format(nodes))


    # Build workflow from nodes
    log.critical('Adding nodes to workflow...')
    wf = parse_nodes(nodes)

    return wf


def _load_data(path):
    with open(path, 'r') as fp:
        return fp.read()


def _load_json(path):
    with open(path, 'r') as fp:
        return json.load(fp)


def _load_yaml(path):
    with open(path, 'r') as fp:
        return yaml.load(fp.read())

