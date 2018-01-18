"""Line-separated text file. Space-delimited 5-column format:
"<chr> <start> <stop> + <target_name>" 1-based coordinates. This format
begins with a sequence dictionary."""

required_fields = ('seqname', 'start', 'stop', 'strand', 'target_name')
optional_fields = None
delimiter = ' '
name = 'Picard-style interval list'
extension = '.interval_list'

#
# def parse_interval(string):
#     """ Given a single interval from a picard interval list file, returns an
#      Interval object"""
#     values = []
#     cols = string.split(delimiter)
#
#     for field in required_fields:
#         values.append((field, cols.pop(0)))
#
#     if cols:
#         # If there are still fields remaining after consuming all
#         # the required and optional fields
#         raise IndexError
#
#     fields = dict(values)
#     return Interval.from_dict(fields)
#
#
#
# def read(path):
#     """This just pops lines from the file and puts them into the sequence
#     dictionary until it reaches a line that does not start with "@" or
#     the end. If there are no sequence_dictionary lines after this process
#     an exception will be raised."""
#     lines = utils.read_lines_allow_gzip(path)
#     interval_file = IntervalFile(format=sys.modules[__name__])
#     sequence_dictionary = []
#
#     while 1:
#         try:
#             line = lines.pop(0)
#             if line.startswith('@'):
#                 sequence_dictionary.append(line)
#             else:
#                 lines.insert(0, line)  # Put it back
#                 break
#         except StopIteration:
#             break
#
#     if not sequence_dictionary:
#         raise Exception('No sequence dictionary found at top of file')
#     else:
#         interval_file.sequence_dictionary = '\n'.join(sequence_dictionary)
#
#     for i, line in enumerate(lines):
#         try:
#             interval = parse_interval(line)
#         except Exception:
#             # Catch any exceptions that occur while parsing the line and
#             # reraise them with line number in the message
#             raise Exception('Reading line {}: {}'.format(i, line))
#         interval_file.intervals.append(interval)
#     return interval_file
