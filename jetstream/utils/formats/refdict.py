"""Picard sequence dictionary format
https://broadinstitute.github.io/picard/command-line-overview.html"""


def read(path):
    with open(path, 'r') as fp:
        data = fp.read()
    return data


def parse_refdict_line(line):
    if line.startswith('@SQ'):
        groups = line.strip('\n').split('\t')[1:]
        groups = [g.split(':', 1) for g in groups]
        fields = {k: v for k, v in groups}
        return fields
    else:
        return None


def refdict_to_bedtools_genome(refdict, out_path):
    """ Converts RefDict format to a Bedtools genome. Returns the
    result as a NamedTemporaryFile """
    refdict = read(refdict)

    with open(out_path, 'w') as fp:
        for line in refdict.splitlines():
            fields = parse_refdict_line(line)
            if fields is not None:
                new_line = "{SN}\t{LN}".format(**fields)
                print(new_line, file=fp)

    return out_path
