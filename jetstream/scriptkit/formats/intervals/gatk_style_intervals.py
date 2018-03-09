"""Line-separated text file. Single column with format <chr>:<start>-<stop>.
 Does not allow for tracks or comments AFAIK. 1-based coordinates:
 https://software.broadinstitute.org/gatk/documentation/article.php?id=1319"""
from .generic import interval

required_fields = ('seqname',)
optional_fields = ('start', 'stop')
metalines = tuple()
delimiter = None
name = 'GATK-style interval list'
extension = '.intervals'


def read(string):
    # TODO sometimes the records are just "chr:start"
    seqname, *others = string.split(':')
    if others:
        start, stop = others.split('-')
    else:
        start = stop = None
    i = interval(seqname=seqname, start=start, stop=stop)
    return i


def write(interval):
    """ Convert an Interval object to a gatk_intervals_list compatible
    interval string """
    if interval['_is_meta']:
        return None
    else:
        string = ''
        seqname = interval.get('seqname')
        if seqname is None:
            msg = 'Required field "{}" is not set.'
            raise ValueError(msg.format(seqname))
        else:
            string += seqname

        start = interval.get('start')
        if start is not None:
            string += ':{}'.format(start)

        stop = interval.get('stop')
        if stop is not None:
            string += '-{}'.format(stop)

        return string
