"""Line separated, tab-delimited, text file. 3 required columns, 9 optional.
Allows intervals to be grouped into tracks by including track-lines. 0-based
coordinates. https://genome.ucsc.edu/FAQ/FAQformat.html"""
import sys

from jetstream import utils
from .intervals import Interval, IntervalFile, Track

required_fields = ('seqname', 'start', 'stop')
optional_fields = ('name', 'score', 'strand', 'thickStart', 'thickEnd',
                   'itemRgb', 'blockCount', 'blockSizes', 'blockStarts')
delimiter = '\t'
name = 'BED (Browser Extensible Data)'
extension = '.bed'
description = __doc__


def parse_interval(string):
    """ Given a single interval from a GFFv2 file, returns an
    Interval object"""
    values = []
    cols = string.split(delimiter)

    for field in required_fields:
        values.append((field, cols.pop(0)))

    try:
        for field in optional_fields:
            values.append((field, cols.pop(0)))
    except IndexError:
        pass

    if cols:
        # If there are still fields remaining after consuming all
        # the required and optional fields
        raise IndexError

    fields = dict(values)
    interval = Interval.from_dict(fields)

    # Account for 0-based indexing
    interval.start += 1
    interval.stop += 1

    return interval


def read(path):
    lines = utils.read_lines_allow_gzip(path)
    interval_file = IntervalFile(format=sys.modules[__name__])
    current_track = None

    for i, line in enumerate(lines):
        try:
            if line.startswith('#'):
                interval_file.comments.append(line)
                continue
            elif line.startswith(('track', 'browser')):
                if interval_file.intervals and current_track is None:
                    # This means that the file has tracks, but there are
                    # intervals outside of the first track. These intervals
                    # are added to a dummy track to preserve their order.
                    blank_track = Track()
                    interval_file.tracks.append(blank_track)
                    blank_track.intervals.append(interval_file.intervals)
                else:
                    if current_track is not None:
                        interval_file.tracks.append(current_track)
                    current_track = Track.from_string(line)
            else:
                interval = parse_interval(line)
                if current_track:
                    current_track.intervals.append(interval)
                interval_file.intervals.append(interval)
        except Exception:
            # Catch any exceptions that occur while parsing the line and
            # reraise them with line number in the message
            raise Exception("Reading line {}: {}".format(i, line))

    return interval_file


def dumps_interval(interval, cols=None):
    res = []
    if cols is not None:
        for field in cols:
            value = getattr(interval, field)
            res.append(value)
    else:
        for field in required_fields:
            value = getattr(interval, field)
            res.append(value)

        for field in optional_fields:
            value = getattr(interval, field, None)
            if value:
                res.append(value)
    return delimiter.join([str(i) for i in res])


def dumps(interval_file, cols=None):
    lines = []
    if interval_file.tracks:
        for track in interval_file.tracks:
            if track:
                lines.append(str(track))
            for interval in track.intervals:
                lines.append(dumps_interval(interval, cols=cols))
    else:
        for interval in interval_file.intervals:
            lines.append(dumps_interval(interval, cols=cols))

    return '\n'.join(lines)
