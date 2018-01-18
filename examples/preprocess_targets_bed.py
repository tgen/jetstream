#!/usr/bin/env python3
import sys
import os
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

    with open(out_path, 'w') as fp:
        print(intervals.to_bed(bed), file=fp)

    return out_path


if __name__ == '__main__':
    bed_path = sys.argv[1]
    out_path = sys.argv[2]

    if os.path.exists(out_path):
        raise OSError("{} already exists".format(out_path))
    else:
        preprocess_bed(bed_path, out_path)



