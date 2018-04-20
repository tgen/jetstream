"""The workflow spec defines the node attributes required/allowed for building
 a workflow graph from a template. """
import logging

log = logging.getLogger(__name__)

SPEC = {
    'name': 'Jetstream node properties definition',
    'spec': {
        'node_properties': {
            'id': {
                'description': 'Unique id for the node',
                'required': True,
                'type': str
            },
            'cmd': {
                'description': 'Array of command arguments',
                'type': list,
                'required': True
            },
            'stdin': {
                'description': 'Data to be sent to the command stdin',
                'type': str
            },
            'stdout': {
                'description': 'Write stdout to a file',
                'type': str
            },
            'stderr': {
                'description': 'Write stderr to a file',
                'type': str
            },
            'stagein': {
                'description': 'Data to be copied into working directory prior '
                               'to execution',
                'type': str
            },
            'stageout': {
                'description': 'Data to be copied out of working directory '
                               '(into final) after execution',
                'type': str
            },
            'before': {
                'description': 'This node should execute before some other '
                               'node/s',
                'type': list,
            },
            'after': {
                'description': 'This node should execute after some other '
                               'node/s',
                'type': list,
            }
        },
    }
}


class NodeValidationError(Exception):
    pass


class SpecNameHere:
    def __init__(self, spec):
        self.spec = spec['spec']
        self.name = spec['name']
        self.node_properties = self.spec['node_properties']
        self.docs = __doc__ # Todo auto format docs with yaml

    def required_props(self):
        for name, attrs in self.node_properties.items():
            if 'required' in attrs and attrs['required'] is True:
                yield name

    def has_unknown_props(self, node):
        errs = []
        for k in node.keys():
            if not k in self.node_properties:
                errs.append('Unknown node property "{}"'.format(k))
        return errs

    def missing_req_props(self, node):
        errs = []
        for prop in self.required_props():
            if not prop in node:
                errs.append('Missing required prop "{}"'.format(prop))
        return errs

    def wrong_prop_type(self, node):
        errs = []
        for k, v in node.items():
            if not k in self.node_properties:
                continue

            cls = self.node_properties[k]['type']
            if not isinstance(v, cls):
                errs.append('Property "{}" should be type "{}"'.format(k, cls))

        return errs

    def validate(self, node):
        errs = []
        errs += self.missing_req_props(node)
        errs += self.has_unknown_props(node)
        errs += self.wrong_prop_type(node)

        if errs:
            msg = '\n'.join(errs)
            log.critical('Node validation failed for {}:\n{}'.format(node, msg))
            return False
        else:
            return True

    def coerce(self, node):
        # Stagein/out can be strings
        if 'stagein' in node:
            if isinstance(node['stagein'], str):
                node['stagein'] = {'src': node['stagein']}

        if 'stageout' in node:
            if isinstance(node['stageout'], str):
                node['stageout'] = {'src': node['stageout']}

        # Any node property that is supposed to be a list can
        # also be a string. The string will be coerced to a
        # len 1 list.
        for k, v in node.items():
            cls = self.node_properties[k]['_class']

            if cls is list and isinstance(v, str):

                log.warning('Coercing attribute "{}: {}" to list.'.format(
                    k, v, cls))
                v = [v, ]
            node[k] = cls(v)


current = SpecNameHere(SPEC)
