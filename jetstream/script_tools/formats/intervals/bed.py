"""Line separated, tab-delimited, text file. 3 required columns, 9 optional.
Allows intervals to be grouped into tracks by including track-lines. 0-based
coordinates. https://genome.ucsc.edu/FAQ/FAQformat.html

Note!: 0-Based Coordinates, every other interval format uses 1-based.
"""
from .generic import interval

required_fields = ('seqname', 'start', 'stop')
optional_fields = ('name', 'score', 'strand', 'thickStart', 'thickEnd',
                   'itemRgb', 'blockCount', 'blockSizes', 'blockStarts')
metalines = ('#', 'track', 'browser')
delimiter = '\t'
name = 'BED (Browser Extensible Data)'
extension = '.bed'
description = __doc__


def read(string):
    """ Given a single interval from a GFFv2 file, returns an Interval object.
    Will return meta lines if they start with #, track, or browser. """
    if string.startswith(metalines):
        return interval(_is_meta=True, seqname=string)

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
        raise IndexError("Too many columns: {}".format(cols))

    fields = dict(values)
    i = interval(**fields)

    # Account for 0-based indexing
    i['start'] += 1
    i['stop'] += 1

    return i


def write(interval):
    interval = interval.copy()
    if interval['_is_meta']:
        if interval['seqname'].startswith(metalines):
            return interval['seqname']
        else:
            return None
    else:
        res = list()
        res.append(interval['seqname'])
        res.append(interval['start'] - 1)  # Account for 0-based coordinates
        res.append(interval['stop'] - 1)

        # TODO: How to handle writing optional fields?
        # for field in optional_fields:
        #     value = interval.get(field)
        #     if value is not None:
        #         res.append(value)

        return delimiter.join([str(i) for i in res])
