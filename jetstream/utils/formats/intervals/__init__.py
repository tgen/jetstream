from jetstream.utils import read_lines_allow_gzip

from . import bed, gffv2, gatk_style_intervals
from .generic import IntervalFile


# TODO test GTF/GFF -> BED dumps
# TODO Do we want to allow for bedtools style functions on generic intervals?
# or should we just convert to bedfile then do it with bedtools?
# TODO GFFv3 IntervalFileFormat
# TODO think about how a coerce(obj, format) function might work when we have
# data that does not meet all the requirements for a format, but want it anyway


def _read_format(path, format):
    lines = [format.read(l) for l in read_lines_allow_gzip(path)]
    return IntervalFile(path, format, lines)


def read_bed(path):
    return _read_format(path, bed)


def read_gffv2(path):
    return _read_format(path, gffv2)


def read_gatk_style_intervals(path):
    return _read_format(path, gatk_style_intervals)


def _to_format(interval_file, format):
    lines = []
    for interval in interval_file:
        line = format.write(interval)
        if line is not None:
            lines.append(line)
    return '\n'.join(lines)


def to_bed(interval_file):
    return _to_format(interval_file, bed)


def to_gffv2(interval_file):
    return _to_format(interval_file, gffv2)


def to_gatk_style_intervals(interval_file):
    return _to_format(interval_file, gatk_style_intervals)

