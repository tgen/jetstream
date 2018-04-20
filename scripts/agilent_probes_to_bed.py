#!/usr/bin/env python
""" Reads the probes file given as first argument and prints a bed file
to stdout. This works for Agilent probe files that:

    - Tab-separated table
    - Includes a "Coordinate" column
    - Coordinates are in the format <chrom>:<start>-<stop>

This also applies some fixes to chromosome names that are commonly malformed
in Agilent files
"""
import re
import sys
import logging

log = logging.getLogger(__name__)

# regex pattern for matching "chrom:start-stop"
pat = re.compile(r'(?P<chrom>[^:]*):(?P<start>[^-]*)-(?P<stop>[^\t]*)')

replace_mistakes = {
    'chrx': 'chrX',
    'chry': 'chrY',
    'chrun_': 'chrUn_'
}

# Read in the lines while using regex to pull useful information
# out of the first column. Prints lines that match the pattern in
# a three-column bed format.
with open(sys.argv[1], 'r') as fp:
    header = fp.readline().strip('\n').split('\t')
    if not 'Coordinates' in header:
        raise ValueError('Unable to find "Coordinates" col in header')

    coords_col_num = header.index('Coordinates')

    for i, line in enumerate(fp):
        try:
            coords_string = line.strip('\n').split('\t')[coords_col_num]
            match = pat.match(coords_string)
        except IndexError:
            log.critical('Not enough fields in line {}: {}'.format(i, line))
            continue

        if match:
            groups = match.groupdict()

            # Fix malformed chromosome names
            for old, new in replace_mistakes.items():
                groups['chrom'] = groups['chrom'].replace(old, new)

            row = "{}\t{}\t{}\t".format(
                groups['chrom'], groups['start'], groups['stop'])

            print(row)

        else:
            log.critical('Unable to match coordinate string pattern for '
                         'line {}: {}'.format(i, line))
