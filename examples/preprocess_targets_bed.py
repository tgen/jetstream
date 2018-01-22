#!/usr/bin/env python3
import sys
import argparse
from jetstream.formats import intervals


def remove_prefix(string, prefix):
    if string.startswith(prefix):
        return string[len(prefix):]
    else:
        return string


def preprocess_bed(path, replace_underscores=False, replace_mt=False,
                   remove_chr=True):
    """ Preprocessing used for canine bed files from Agilent """
    bed = intervals.read_bed(path)

    # Iterate over every interval in the bed file
    for i in bed:
        if remove_chr:
            i['seqname'] = remove_prefix(i['seqname'], 'chr')
        if replace_underscores:
            i['seqname'] = i['seqname'].replace('_', '.')
        if replace_mt:
            if i['seqname'] == 'M':
                i['seqname'] = 'MT'

    print(intervals.to_bed(bed))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Preprocess bed files prior to exome kit refpack creation. "
                    "This will remove 'chr' prefixes from chromosome names and "
                    "replace chromosome 'M' with 'MT'"
    )

    parser.add_argument('bed', help='path to a bed file to preprocess')

    parser.add_argument('--replace-underscores', action='store_true', default=False,
                        help="Replace underscores with dots. Useful for some "
                             "Agilent files that don't match standard chrom "
                             "names, but can cause issues on most others.")

    parser.add_argument('--replace-mt', action='store_true', default=False,
                        help="Replace chromosome name 'M' with 'MT' to match "
                             "reference genomes.")

    parser.add_argument('--dont-remove-chr',
                        dest='remove_chr', action='store_false', default=True,
                        help="Don't remove the chr prefix from chromosome "
                             "names")

    args = parser.parse_args()

    preprocess_bed(
        path=args.bed,
        replace_underscores=args.replace_underscores,
        replace_mt=args.replace_mt,
        remove_chr=args.remove_chr,
    )



