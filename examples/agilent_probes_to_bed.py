#!/usr/bin/env python
""" Reads the probes file given as first argument and prints a bed file to stdout """
import re
import sys

# regex pattern for matching "chrom:start-stop"
pat = re.compile(r'(?P<chrom>[^:]*):(?P<start>[^-]*)-(?P<stop>[^\t]*)')

# Read in the lines while using regex to pull useful information
# out of the first column. Prints lines that match the pattern in
# a three-column bed format.
with open(sys.argv[1], 'r') as fp:
    header = fp.readline()
    for line in fp:
        match = pat.match(line.split('\t')[0])
        if match:
            groups = match.groupdict()
            if groups['chrom'] == 'chrx':
                groups['chrom'] = 'chrX'
            row = "{}\t{}\t{}\t".format(groups['chrom'], groups['start'], groups['stop'])
            print(row)
