""" Generalized representation of intervals that can be converted to different
formats. Interval formats need to define read and write functions.

Read should take a string as its only argument and return an interval dict or
raise an index error if unable to parse the interval string.

Write should take an interval dict as its only argument and return an interval
string or raise an index error if the interval lacks necessary fields. The
writer may return None in order to skip metadata lines.

Both methods should recognize the "_is_meta" field and handle those intervals
appropriately. Some formats allow metadata: comments, tracks lines, or browser
lines. If a format allows these metadata, the reader should accept them and the
writer should return them.

"""


def interval(_is_meta=False, **kwargs):
    """ All intervals require at least a seqname, start, and stop. Other
    properties are optional and can be given as kwargs. """
    if _is_meta:
        kwargs['_is_meta'] = True
        kwargs['seqname'] = str(kwargs['seqname'])
    else:
        kwargs['_is_meta'] = False
        kwargs['seqname'] = str(kwargs['seqname'])
        kwargs['start'] = int(kwargs['start'])
        kwargs['stop'] = int(kwargs['stop'])
    return kwargs


class IntervalFile(object):
    def __init__(self, path, format, intervals):
        self.format = format
        self.path = path
        self.intervals = intervals

    def __getitem__(self, item):
        return self.intervals.__getitem__(item)

    def __getattr__(self, item):
        return self.intervals.__getattribute__(item)

    def __len__(self):
        return len(self.intervals)

    def filter(self, fn, meta=False):
        """ Apply a filter fn to all intervals in this file, returns
        a new IntervalFile. This always includes meta intervals unless
        "meta" is set, then the filter fn is also applied to those
        intervals """
        if meta:
            intervals = [i.copy() for i in self.intervals if fn(i)]
        else:
            intervals = []
            for i in self.intervals.copy():
                if i['_is_meta']:
                    intervals.append(i.copy())
                elif fn(i):
                    intervals.append(i.copy())

        return IntervalFile(self.path, self.format, intervals)

    def map(self, fn, meta=False):
        """ Map a function to all intervals in this file, returns a
        list of results. Meta intervals will be skipped unless "meta"
        is True."""
        if meta:
            results = [fn(i) for i in self.intervals.copy()]
        else:
            results = []
            for i in self.intervals:
                if i['_is_meta']:
                    pass
                else:
                    results.append(fn(i))
        return results
