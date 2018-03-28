"""
# Portability

This currently only works on Linux machines locally or with a Slurm
batch scheduler. It also relies heavily on Python.

# Upfront workflow rendering

Workflows are an immutable data structure. Rendering a workflow
is a dynamic process that responds to input data and template
directives. But, the workflow that is generated should be a
complete description of the tasks required, and the order of
execution.

The advantages of this are that complete workflows can be built
prior to runtime. You can export exact diagrams, directory structures,
command lists, etc.. without executing any steps.

What about feedback? Conditionals?

The workflow is an immutable network graph defined prior to runtime.
It cannot be modified by events that occur during runtime. Data can
be captured from nodes, and stored in the project data. Downstream
nodes can then make use of the project data for evaluating expression
attributes. This is the only feedback mechanism.

If a node exists in a workflow, the runner will always launch it. The
only exceptions to this rule occur if the runner exits before reaching
that node, or one of its dependencies has failed. Nodes can fail
gracefully by condtional expressions. In other words, nodes can be
skipped in some cases. But the runner will still lauch the node, and
it will still fail, but it will not count towards a workflow failure.

Cases where feedback may be necessary:

Dynamic chunking of data

Some input data needs to be split into n chunks where n is determined
during workflow runtime. Each chunk then needs to be treated as an
individual task in the workflow by downstream tasks.

Note that this is not a problem if n can be determined prior to runtime,
or if your application can handle the chunking internally.


# Parsing/Serializing

json, yaml, etc..

"""
import logging

log = logging.getLogger(__name__)
SPEC = {
    'name': '<spec name here> Definition',
    'version': '0.1a',
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
                'description': 'Capture stdout to a project variable',
                'type': str
            },
            'stderr': {
                'description': 'Capture stderr to project variable',
                'type': str
            },
            'stagein': {
                'description': 'Data to be copied into working directory prior '
                               'to execution',
                'type': str
            },
            'stageout': {
                'description': 'Data to be copied out of workfing directory '
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
        self.version = spec['version']
        self.node_properties = self.spec['node_properties']
        self.docs = __doc__ # Todo auto format docs with yaml
        log.critical('Spec {} loaded'.format(self.version))


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
