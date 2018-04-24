"""Workflow builder module contains functions required for building a workflow
object from a workflow template (yaml file)"""
import logging
from jetstream import utils
from jetstream.workflows import Workflow
from jetstream.workflows import spec

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


def build_workflow(nodes):
    """Given a sequence of nodes (dictionaries with properties described in
     the workflow specification), returns a workflow with nodes and edges
     built. """
    wf = Workflow()

    for node in nodes:
        log.debug('Adding node to workflow: {}'.format(node))
        wf.add_node(node_id=node['id'], **node)

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
