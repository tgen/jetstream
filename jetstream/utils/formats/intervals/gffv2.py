"""Line-separated, tab-delimited, text file. 1-based coordinates, 9 required
columns, 0 optional columns. Allows intervals to be grouped into tracks by
including track-lines. https://www.ensembl.org/info/website/upload/gff.html"""
from .generic import interval

required_fields = ('seqname', 'source', 'feature', 'start', 'stop', 'score',
                   'strand', 'frame', 'attribute')
optional_fields = tuple()
metalines = ('#', 'track')
delimiter = '\t'
name = 'GFF (General Feature Format) v2'
extension = '.gff'
description = __doc__


def read(string):
    """ Converts a single GFFv2 compatible interval string to an Interval
    object"""
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
        raise IndexError

    fields = dict(values)
    return interval(**fields)


def write(interval):
    """ Convert an Interval object to a GFFv2 compatible interval string """
    if interval['_is_meta']:
        if interval['seqname'].startswith(metalines):
            return interval['seqname']
        else:
            return None
    else:
        cols = []
        for field in required_fields:
            value = interval.get(field)
            if value is not None:
                cols.append(value)
            else:
                msg = 'Required field "{}" value is not set.'
                raise ValueError(msg.format(field))

        for field in optional_fields:
            value = interval.get(field)
            if value is not None:
                cols.append(value)

        return delimiter.join([str(i) for i in cols])
