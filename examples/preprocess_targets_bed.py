#!/usr/bin/env python3
import sys
from jetstream.formats import intervals


def remove_prefix(string, prefix):
    if string.startswith(prefix):
        return string[len(prefix):]
    else:
        return string


def preprocess_bed(path, out_path):
    """ Preprocessing used for canine """
    bed = intervals.read_bed(path)

    # Remove chr prefix from sequence names
    for i in bed:
        i['seqname'] = remove_prefix(i['seqname'], 'chr')
        i['seqname'] = i['seqname'].replace('_', '.')

    print(intervals.to_bed(bed))


if __name__ == '__main__':
    bed_path = sys.argv[1]
    preprocess_bed(bed_path)



