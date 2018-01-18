#!/usr/bin/env python3
import sys
import tempfile
from jetstream.formats import intervals

ints = intervals.read_gffv2(sys.argv[1])
ints = ints.filter(lambda i: i['feature'] == 'exon')  # Select only exons
exon_bed = tempfile.NamedTemporaryFile()
with open(exon_bed.name, 'w') as fp:
    fp.write(intervals.to_bed(ints))

