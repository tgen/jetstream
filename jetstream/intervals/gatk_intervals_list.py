"""Line-separated text file. Single column with format <chr>:<start>-<stop>.
 Does not allow for tracks or comments AFAIK. 1-based coordinates:
 https://software.broadinstitute.org/gatk/documentation/article.php?id=1319"""
import sys

from jetstream import utils
from .intervals import Interval, IntervalFile

required_fields = ('seqname',)
optional_fields = ('start', 'stop')
delimiter = None
name = 'GATK-style interval list'
extension = '.intervals'


def parse_interval(string):
    seqname, *others = string.split(':')
    if others:
        start, stop = others.split('-')
    else:
        start = stop = None
    interval = Interval(seqname, start, stop)
    return interval


def read(path):
    lines = utils.read_lines_allow_gzip(path)
    interval_file = IntervalFile(format=sys.modules[__name__])

    for i, line in enumerate(lines):
        try:
            interval = parse_interval(line)
        except Exception:
            # Catch any exceptions that occur while parsing the line and
            # reraise them with line number in the message
            raise Exception('Reading line {}: {}'.format(i, line))
        interval_file.intervals.append(interval)

    return interval_file
