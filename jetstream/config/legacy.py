"""
# Comments from original validation script:

### check if we can process this config or just hold it
### validate email address
### validate number of fastq lines > 0
### validate fastq files exist (R2 also)
### validate fastq read group exists
### validate project name
### validate results directory is present
### validate pipeline name is present
# Pipeline name is valid
### validate recipe name
# Recipe name is valid
### validate save recipe
# save recipe is valid
### validate debit account name
#validate sample is only defined once
#validate sample name does not contain a period
#validate kit code is 5 characters
#validate kit code for exomes for this sample is valid
### validate kit name, sample name
# validate RNA strandedness is known
# validate RNA has star index defined for this read length
#validate triplet line contains sample names that are valid

### check if genome samples exist
    ### if allele count is requested make sure first two is set up as a dna pair
    ### and make sure third variable is RNA
    # check  that CIRCOS variable is set
    # if circos is yes we have to make sure circos exome pair exists and theyre
    both exome
### validate if RNA sample exists TOPHATGTF, USEMASK vars are defined

### validate FQ= in Pegasus pipe formatted correctly

### Validate dnaPair samples exist
### Validate dnaFami samples exist
### Validate rnaPair samples exist
"""
import logging
from collections import OrderedDict

from jetstream.utils import Source

log = logging.getLogger(__name__)


class ConfigParsingException(Exception):
    """ Raised when config files contain syntax errors """


def deserialize(data):
    """Loads a legacy config file from data. Performs  several validation
    steps, briefly:

     - Leading and trailing whitespace is removed
     - Checks for single occurrence of "=START" and "=END" lines
     - Verifies that no data follows "=END"
     - All metadata keys are translated to lowercase

    """
    source = Source(data.strip())
    meta_lines, sample_lines = _split_sections(source)
    meta = _parse_meta_lines(meta_lines)
    sample_line_groups = _group_sample_lines(sample_lines)
    samples = [_parse_sample(group) for group in sample_line_groups]
    return {'meta': meta, 'samples': samples}


def read(path):
    with open(path, 'r') as fp:
        return deserialize(fp.read())


def _split_sections(source):
    """ Split the metadata lines from the sample lines """
    lines = source.splitlines()

    start_str = '=START'
    end_str = '=END'

    if lines.count(start_str) != 1:
        raise ConfigParsingException('Error finding "=START" line')
    else:
        start = lines.index('=START')

    if lines.count(end_str) != 1:
        raise ConfigParsingException('Error finding "=END" line')
    else:
        end = lines.index('=END')

    if lines[-1] != "=END":
        raise ConfigParsingException('Lines found after "=END"')

    meta = lines[0: start]
    sample = lines[start + 1: end]
    return meta, sample


def _parse_meta_lines(lines):
    """ Each line in metadata should follow a "key=value" syntax. This
    function splits the lines in metadata into key-value pairs """
    kv_tuples = []
    for line in lines:
        try:
            key, value = line.split('=')
            key = key.lower()
            kv_tuples.append((key, value))
        except ValueError:
            msg = 'Error parsing metadata: "{}"'.format(line.print_ln())
            raise ConfigParsingException(msg)
    meta = OrderedDict(kv_tuples)
    return meta


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
    """Create a structured sample from a list of sample lines

    This diverges from the old structure a little bit:

    The kit code, is associated with a sample in the config files. But the
    assay is not really a property of the sample itself, it's a property of
    the data object. For instance, one sample can be prepped multiple times
    with different kits. So, here kit becomes a property of each data object
    """
    sample_line = lines[0]
    sample = {}
    sample['data'] = []
    sample['line_number'] = sample_line.line_number

    # Get info from the sample line first
    fields = sample_line.strip('SAMPLE=').split(',')  # Break it up
    kit = fields[0]  # this becomes a property of the data objs
    sample['name'] = fields[1]
    assay = fields[2]
    if len(fields) > 3:
        sample['dilution_id'] = fields[3]

    # Then parse each data object listed for that sample
    for line in lines[1:]:
        data = {}
        data['line_number'] = line.line_number
        data['type'], fields = line.split('=') # ex: "FQ=..." type is FQ

        # Split the remaining fields on commas
        fields = fields.split(',')
        data['read_group'] = fields[0]
        data['path'] = fields[1]
        # TODO Sometimes there's a third or fourth field I think?
        data['kit'] = kit
        data['assay'] = assay
        sample['data'].append(data)

    return sample
