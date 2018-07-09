"""Tools for parsing legacy config files and converting to new formats."""
import os
import logging

from jetstream.utils import Source, records_to_csv, yaml

log = logging.getLogger(__name__)


class ConfigParsingException(Exception):
    """ Raised when config files contain syntax errors """


def load(path):
    """Reads a legacy config file from the given path, returns
    config object"""
    with open(path, 'r') as fp:
        return loads(fp.read())


# def dump(*args, **kwargs):
#     """ There is no dumper for legacy config files. You shouldn't be
#      converting to the old format """
#     raise NotImplementedError


def loads(data):
    """Generates a legacy config file from data. Performs several validation
    steps and adds some additional attributes. Briefly:

     - Leading and trailing whitespace is removed
     - Pipeline name is either 'medusa' or 'pegasus' (case insensitive)
     - Contains single occurrence of "=START" and "=END" lines
     - No data following "=END" line
     - All run_parameters keys are translated to lower-case
     - "Kit" given in a SAMPLE line is also added to each data object for
       that sample
     - "Assay" value given in a SAMPLE line is also added to each data object
       for that sample

    """
    source = Source(data.strip())

    # Everything before =START is considered the run_parameters lines,
    # everything after is considered the sample lines.
    run_parameters_lines, sample_lines = _split_sections(source)

    run_parameters = _parse_run_parameters_lines(run_parameters_lines)
    which_pipeline = run_parameters['pipeline'].casefold()
    if which_pipeline not in ('pegasus', 'medusa'):
        raise ConfigParsingException(
            'Unknown pipeline: {}'.format(run_parameters['pipeline']))

    # Parse sample lines
    sample_line_groups = _group_sample_lines(sample_lines)
    data_objects = []
    sample_objects = []
    for group in sample_line_groups:
        try:
            datas, sample = _parse_sample(group)

            # datas is a list of data objects we want to extend
            # data_objects to include
            data_objects.extend(datas)

            # sample there is only one sample so just append to
            # sample_objects
            sample_objects.append(sample)

        except Exception as err:
            msg = '{} while parsing group:\n{}'.format(err, '\n'.join(group))
            raise ConfigParsingException(msg)

    # Add rg and assumed read style
    data_objects = _add_extra_properties(data_objects)

    return {
        'run_parameters': run_parameters,
        'data': data_objects,
        'samples': sample_objects
    }


def explode(config, outdir=''):
    if outdir:
        os.mkdir(outdir)

    samples_path = os.path.join(outdir, 'samples.csv')
    data_path = os.path.join(outdir, 'data.csv')
    run_parameters_path = os.path.join(outdir, 'run_parameters.yaml')

    records_to_csv(config['samples'], samples_path)
    records_to_csv(config['data'], data_path)

    if os.path.exists(run_parameters_path):
        raise FileExistsError(run_parameters_path)
    with open(run_parameters_path, 'w') as fp:
        yaml.dump(config['run_parameters'], stream=fp)


def _split_sections(source):
    """ Split the run_parameters lines from the sample lines """
    lines = [line for line in source.splitlines() if line]

    start_line = '=START'
    end_line = '=END'

    # Check that start_str exists once and only once
    if lines.count(start_line) != 1:
        raise ConfigParsingException('No "=START" line found!')
    else:
        start = lines.index('=START')

    # Check that end_str exists once and only once
    if lines.count(end_line) != 1:
        raise ConfigParsingException('No "=END" line found!')
    else:
        end = lines.index('=END')

    # Check that nothing is present after end line
    if lines[-1] != "=END":
        raise ConfigParsingException('Lines found after "=END"')

    # Now split the lines into two groups, run_parameters and sample
    run_parameters = lines[0: start]
    sample = lines[start + 1: end]
    return run_parameters, sample


def _parse_run_parameters_lines(lines):
    """ Each line in run_parameters should follow a "key=value" syntax. This
    function splits the lines in run_parameters into key-value pairs and
    returns an OrderedDict of the run_parameters.

    Note: This ignores case on the keys and all keys in the resulting
    dict will be converted to lower-case.

    """
    run_parameters = dict()

    kv_tuples = list()
    for line in lines:
        try:
            key, value = line.split('=')
            kv_tuples.append((key.lower(), value))
        except ValueError:
            msg = 'Error parsing run_parameters: "{}"'.format(line.print_ln())
            raise ConfigParsingException(msg) from None

    # NOTE: Here I could have dynamically set value types to
    # lists if the key was found multiple times. But, This could have
    # potentially masked user errors in the config files. I've chosen
    # to hard-code the fields which can be lists, in order to make the
    # behavior more predictable. Every other value is parsed with str()

    for k, v in kv_tuples:
        if k in ('jirset', 'dnapair', 'triplet4allelecount'):
            if k not in run_parameters:
                run_parameters[k] = list()
            run_parameters[k].append(v)
        else:
            run_parameters[k] = str(v)

    return run_parameters


# The following functions convert the sample lines in a config file into
# data objects (nested dictionary/list structure). This causes a couple small
# changes in the structure:
#
#    1) Data objects inherit 'kit' and 'assay' attributes from their SAMPLE
#
#    The kit code and assay are associated with a sample in the config files.
#    But, the assay is not really a property of the sample itself, it's a
#    property of the data object. For instance, one sample can be prepped
#    multiple times with different kits. Here, kit is copied as a property of
#    each data object.
#
#
#    2) Data will ALWAYS have a library property.
#
#    Dilution ID is part of the sample data for Medusa config files. This
#    was later changed in Pegasus. The dilution ID is also called a library,
#    and is eventually used for making read groups. This parser always adds
#    dilution id (library) as a property of the data.
#
def _group_sample_lines(lines):
    """ Groups sample lines together with data lines """
    samples = []

    sample_iter = iter(lines)
    line = next(sample_iter)

    if line.startswith('SAMPLE='):
        current_group = [line]
    else:
        raise ConfigParsingException('Sample lines should start with SAMPLE')

    while 1:
        try:
            line = next(sample_iter)
        except StopIteration:
            if current_group:
                samples.append(current_group)
            break

        if line.startswith('SAMPLE='):
            samples.append(current_group)
            current_group = [line]
        else:
            current_group.append(line)

    return samples


def _parse_sample(lines):
    """Create a list of data objects from a list of sample lines"""
    data_objects = list()
    sample_line = lines[0]
    s = {'line_number': sample_line.line_number}

    # The sample fields are everything after 'SAMPLE=' Collect data
    # from the sample line and add it to the sample dict.
    fields = sample_line.partition('SAMPLE=')[2]
    try:
        # Pegasus field order
        s['kit'], s['sample_name'], s['assay'], s['library'] = fields.split(',')
    except ValueError:
        # Medusa field order
        s['kit'], s['sample_name'], s['assay'] = fields.split(',')

    # Every line after 'SAMPLE=' is a data object
    # Collect data information from the data line in a dictionary
    # and append it to sample['data'] array
    for line in lines[1:]:
        d = {
            'line_number': line.line_number,
            'sample_name': s['sample_name'],
            'kit': s['kit'],
            'assay': s['assay']
        }
        d['type'], _, data_fields = line.partition('=')  # ex: "FQ=..."
        d['rg_id'], d['path'] = data_fields.split(',')

        try:
            # Pegasus has 3 fields in the rg field
            d['fcid'], d['lane'], d['library'] = d['rg_id'].split('_')
        except ValueError:
            # Medusa has only 2
            d['fcid'], _, d['lane'] = d['rg_id'].partition('_')
            d['library'] = s['library']

        data_objects.append(d)

    return data_objects, s


def _add_extra_properties(data):
    """Copies some data from the sample properties to all of its data
    objects. Also generates read group data. Note this follows the TGen
    convention of PU: FCID_LANE and ID:FCID_LANE[_LIBRARY] where LIBRARY
    may not be present in all RG:ID tags. This is generally unimportant
    unless you are considering merging SAMs. """
    for obj in data:
        if obj['type'] == 'FQ':
            obj['read_style'] = 'paired-end'

        obj['read_group'] = {
            'ID': obj['rg_id'],
            'CN': 'TGen',
            'LB': obj['library'],
            'PL': 'ILLUMINA',
            'PU': '{}_{}'.format(obj['fcid'], obj['lane']),
            'SM': obj['sample_name']
        }

    return data
