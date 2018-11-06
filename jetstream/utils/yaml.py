"""Monkey patch the default yaml dump to prefer block style for sequences

Notes: PyYaml is currently undergoing some extensive changes aimed to
improve the API. They are hoping to release along with 3.7. It may provide
a mechanism for setting default_flow_style at a module level.
But, this works for now.

"""
import yaml


def represent_none(dumper, _):
    """Configures the yaml engine to represent null values as blanks"""
    return dumper.represent_scalar('tag:yaml.org,2002:null', '')


def represent_str(dumper, data):
    if len(data.splitlines()) > 1:  # check for multiline string
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
    return dumper.represent_scalar('tag:yaml.org,2002:str', data)


yaml.add_representer(str, representer=represent_str)
yaml.add_representer(type(None), represent_none)
