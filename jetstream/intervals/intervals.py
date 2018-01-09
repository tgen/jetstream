import sys
import gzip

# TODO test GTF/GFF -> BED dumps
# TODO Allow merging/concat interval methods
# TODO create dumps(obj, format) generic function
# TODO create dump(obj, format, fp) generic functions and methods
# TODO GFFv3 IntervalFileFormat
# TODO think about how a coerce(obj, format) function might work


class Interval(object):
    """ Intervals require at least a seqname, start, and stop. Other
    properties are optional and can be given as kwargs. Some recognized 
    properties are already set to None """
    def __init__(self, seqname, start, stop, source=None, feature=None,
                 score=None, strand=None, frame=None, attribute=None,
                 thickStart=None, thickEnd=None, itemRgb=None, blockCount=None,
                 blockSizes=None, blockStarts=None, **kwargs):
        self.seqname = str(seqname)
        self.start = int(start)
        self.stop = int(stop)
        self.source = source
        self.feature = feature
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attribute = attribute
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts

        for k, v in kwargs.items():
            setattr(self, k, v)

    def __str__(self):
        return str(self.__dict__)

    def __len__(self):
        return self.stop - self.start

    @staticmethod
    def from_dict(d):
        return Interval(**d)

    def pad(self, n):
        self.start = max((self.start - n), 1)  # Don't allow start less than 1
        self.stop = self.stop + n


class Track(object):
    def __init__(self, *args):
        """ Properties given as args will be stored in Track.args. If no
        properties are delcared Track.args will be None. Tracks with args
        set to None evaluate to False and will not be included as a line in
        interval files when dumping to string/files. This allows dummy tracks
        to be created to preserve the order of groups of intervals in an
        IntervalFile object. """
        self.intervals = []

        if args:
            self.args = args
        else:
            self.args = None

    def __str__(self):
        if self.args is None:
            return ''
        else:
            return ' '.join(self.args)

    def __nonzero__(self):
        return bool(self.args)

    @staticmethod
    def from_string(string):
        """Create an Track object from a track string. """
        if not string.startswith(('browser', 'track')):
            raise ValueError('"{}" doesnt start with "browser" or "track"'.format(string))
        track = Track(string)
        return track


class IntervalFile(object):
    def __init__(self, format, path=None, tracks=None, intervals=None, comments=None):
        self.format = format
        self.path = path
        self.tracks = tracks or []
        self.intervals = intervals or []
        self.comments = comments or []
        self.sequence_dictionary = None

#
# def read(path, format):
#     """ Read an interval file. Format should be a IntervalFileFormat """
#     if isinstance(format, IntervalFileFormat):
#         return format.read(path)
#     elif isinstance(format, str):
#         cls = _lookup_class(format)
#         return cls.read(path)
#
#
# def dumps(obj, format):
#     """ Dump an IntervalFile object to a string. Format should be an
#     IntervalFileFormat """
#     if isinstance(format, IntervalFileFormat):
#         return format.dumps(obj)
#     elif isinstance(format, str):
#         cls = _lookup_class(format)
#         return cls.dumps(obj)
#
#
# def _lookup_class(name):
#     return getattr(sys.modules[__name__], name)
#


# def read(path, file_format):
#     if file_format is None:
#         file_format = guess_file_format(path)
#     return _load_interval_file(path, file_format)


# def guess_file_format(path):
#     # TODO
#     return Bed



# def coerce_sequence_names(bed, seqnames):
#     newbed = []
#     for line in bed.splitlines():
#         seqname, start, end = line.split('\t')
#         seqname = seqname.replace('chr',)
#         newbed.append()


# def make_bed(reference_gtf, if_feature='exon'):
#     """ Generate a bed file containing extracted feature coordinates from
#     the reference_gtf, returns bed file as string """
#     bed = []
#     if is_gzip(reference_gtf):
#         with gzip.open(reference_gtf, 'rb') as fp:
#             lines = [line.decode('utf8') for line in fp.readlines()]
#     else:
#         with open(reference_gtf, 'r') as fp:
#             lines = fp.readlines()

#     for i, line in enumerate(lines):
#         line = line.strip()
#         if line == '' or line.startswith('#'):
#             continue
#         else:
#             try:
#                 seqname, source, feature, start, end, score, strand, frame, attribute = line.split('\t')
#             except ValueError as e:
#                 raise Exception('unpacking line {}: {}'.format(i, line))

#         if feature == if_feature:
#             bed.append('\t'.join((seqname, start, end)))
#         else:
#             continue

#     bed = '\n'.join(bed + ['\n'])
#     return bed

