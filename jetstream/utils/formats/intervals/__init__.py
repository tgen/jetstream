"""Intervals are a generalization of several file formats:

GFF
====

All GFF formats (GFF2, GFF3 and GTF) are tabular files with 9 fields per line,
separated by tabs. They all share the same structure for the first 7 fields,
while differing in the definition of the eighth field and in the content and
format of the ninth field. The general structure is as follows:


GFFv3
------

Amazing description of the format:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

"Although there are many richer ways of representing genomic features via XML
and in relational database schemas, the stubborn persistence of a variety of
ad-hoc tab-delimited flat file formats declares the bioinformatics community's
need for a simple format that can be modified with a text editor and processed
with shell tools like grep. "


http://gmod.org/wiki/GFF3
~~~~~~~~~~~~~~~~~~~~~~~~~~

GFF3 Format
GFF3 format is a flat tab-delimited file. The first line of the file is a
comment that identifies the file format and version. This is followed by a
series of data lines, each one of which corresponds to an annotation.Here is a
miniature GFF3 file:

::

    ##gff-version 3
    ctg123  .  exon  1300  1500  .  +  .  ID=exon00001
    ctg123  .  exon  1050  1500  .  +  .  ID=exon00002
    ctg123  .  exon  3000  3902  .  +  .  ID=exon00003
    ctg123  .  exon  5000  5500  .  +  .  ID=exon00004
    ctg123  .  exon  7000  9000  .  +  .  ID=exon00005

The ``##gff-version 3`` line is required and must be the first line of the
file. It introduces the annotation section of the file.

The 9 columns of the annotation section are as follows:

Column 1: "seqid"

The ID of the landmark used to establish the coordinate system for the current
feature. IDs may contain any characters, but must escape any characters not in
the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not contain unescaped
whitespace and must not begin with an unescaped ">".
To escape a character in this, or any of the other GFF3 fields, replace it with
the percent sign followed by its hexadecimal representation. For example, ">"
becomes "%E3". See URL Encoding (or: 'What are those "%20" codes in URLs?') for
details.

Column 2: "source"

The source is a free text qualifier intended to describe the algorithm or
operating procedure that generated this feature. Typically this is the name of
a piece of software, such as "Genescan" or a database name, such as "Genbank."
In effect, the source is used to extend the feature ontology by adding a
qualifier to the type creating a new composite type that is a subclass of the
type in the type column. It is not necessary to specify a source. If there is
no source, put a "." (a period) in this field.

Column 3: "type"

The type of the feature (previously called the "method"). This is constrained
to be either: (a) a term from the "lite" sequence ontology, SOFA; or (b) a SOFA
accession number. The latter alternative is distinguished using the syntax
SO:000000. This field is required.

Columns 4 & 5: "start" and "end"

The start and end of the feature, in 1-based integer coordinates, relative to
the landmark given in column 1. Start is always less than or equal to end.
For zero-length features, such as insertion sites, start equals end and the
implied site is to the right of the indicated base in the direction of the
landmark. These fields are required.

Column 6: "score"

The score of the feature, a floating point number. As in earlier versions of
the format, the semantics of the score are ill-defined. It is strongly
recommended that E-values be used for sequence similarity features, and that
P-values be used for ab initio gene prediction features. If there is no score,
put a "." (a period) in this field.

Column 7: "strand"

The strand of the feature. + for positive strand (relative to the landmark), -
for minus strand, and . for features that are not stranded. In addition, ? can
be used for features whose strandedness is relevant, but unknown.

Column 8: "phase"

For features of type "CDS", the phase indicates where the feature begins with
reference to the reading frame. The phase is one of the integers 0, 1, or 2,
indicating the number of bases that should be removed from the beginning of this
feature to reach the first base of the next codon. In other words, a phase of
"0" indicates that the next codon begins at the first base of the region
described by the current line, a phase of "1" indicates that the next codon
begins at the second base of this region, and a phase of "2" indicates that
the codon begins at the third base of this region. This is NOT to be confused
with the frame, which is simply start modulo 3. If there is no phase, put a
"." (a period) in this field.

For forward strand features, phase is counted from the start field. For reverse
strand features, phase is counted from the end field.
The phase is required for all CDS features.

Column 9: "attributes"

A list of feature attributes in the format tag=value. Multiple tag=value pairs
are separated by semicolons. URL escaping rules are used for tags or values
containing the following characters: ",=;". Spaces are allowed in this field,
but tabs must be replaced with the %09 URL escape. This field is not required.

Column 9 Tags

Column 9 tags have predefined meanings:

ID
Indicates the unique identifier of the feature. IDs must be unique within the
scope of the GFF file.

Name
Display name for the feature. This is the name to be displayed to the user.
Unlike IDs, there is no requirement that the Name be unique within the file.

Alias
A secondary name for the feature. It is suggested that this tag be used
whenever a secondary identifier for the feature is needed, such as locus names
and accession numbers. Unlike ID, there is no requirement that Alias be unique
within the file.

Parent
Indicates the parent of the feature. A parent ID can be used to group exons
into transcripts, transcripts into genes, and so forth. A feature may have
multiple parents. Parent can *only* be used to indicate a partof relationship.

Target
Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide
alignment. The format of the value is "target_id start end [strand]", where
strand is optional and may be "+" or "-". If the target_id contains spaces,
they must be escaped as hex escape %20.

Gap
The alignment of the feature to the target if the two are not collinear
(e.g. contain gaps). The alignment format is taken from the CIGAR format
described in the Exonerate documentation.
http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate?cvsroot=Ensembl).
See the GFF3 specification for more information.

Derives_from
Used to disambiguate the relationship between one feature and another when the
relationship is a temporal one rather than a purely structural "part of" one.
This is needed for polycistronic genes. See the GFF3 specification for more
information.

Note
A free text note.

Dbxref
A database cross reference. See the GFF3 specification for more information.

Ontology_term
A cross reference to an ontology term. See the GFF3 specification for more
information.

Multiple attributes of the same type are indicated by separating the values
with the comma "," character, as in:

Parent=AF2312,AB2812,abc-3

Note that attribute names are case sensitive. "Parent" is not the same as
"parent".

All attributes that begin with an uppercase letter are reserved for later use.
Attributes that begin with a lowercase letter can be used freely by
applications. You can stash any semi-structured data into the database by using
one or more unreserved (lowercase) tags.


GTF
----

The Gene transfer format (GTF) is a refinement of GFF Version 2 and is
sometimes referred to as GFF2.5.[1]

http://genome.ucsc.edu/FAQ/FAQformat.html#format4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GTF (Gene Transfer Format, GTF2.2) is an extension to, and backward compatible
with, GFF2. The first eight GTF fields are the same as GFF. The feature field
is the same as GFF, with the exception that it also includes the following
optional values: 5UTR, 3UTR, inter, inter_CNS, and intron_CNS. The group field
has been expanded into a list of attributes. Each attribute consists of a
type/value pair. Attributes must end in a semi-colon, and be separated from
any following attribute by exactly one space.
The attribute list must begin with the two mandatory attributes:
gene_id value - A globally unique identifier for the genomic source of the
sequence.
transcript_id value - A globally unique identifier for the predicted transcript.

Example:
Here is an example of the ninth field in a GTF data line:

::

    gene_id "Em:U62317.C22.6.mRNA"; transcript_id "Em:U62317.C22.6.mRNA"; exon_number 1

The Genome Browser groups together GTF lines that have the same transcript_id
value. It only looks at features of type exon and CDS.
For more information regarding the GTF2.2 UCSC supported format, see
http://mblab.wustl.edu/GTF22.html. If you would like to obtain browser data in
GTF format, please refer to Genes in gtf or gff format on the wiki.


GFFv2
-----

There are several differing descriptions of the GFFv2 spec. Most seem to
reference GMOD (http://gmod.org/wiki/GFF2) as the authoritative source:


http://gmod.org/wiki/GFF2
~~~~~~~~~~~~~~~~~~~~~~~~~~

The GFF2 File Format
The GFF format is a flat tab-delimited file, each line of which corresponds to
an annotation, or feature. Each line has nine columns and looks like this:

::

    Chr1  curated  CDS 365647  365963  .  +  1  Transcript "R119.7"

The 9 columns are as follows:

reference sequence

This is the ID of the sequence that is used to establish the coordinate system
of the annotation. In the example above, the reference sequence is "Chr1".

source

The source of the annotation. This field describes how the annotation was
derived. In the example above, the source is "curated" to indicate that the
feature is the result of human curation. The names and versions of software
programs are often used for the source field, as in "tRNAScan-SE/1.2".

method

The annotation method, also known as type. This field describes the type of the
annotation, such as "CDS". Together the method and source describe the
annotation type.

start position

The start of the annotation relative to the reference sequence.

stop position

The stop of the annotation relative to the reference sequence. Start is always
less than or equal to stop.

score

For annotations that are associated with a numeric score (for example, a
sequence similarity), this field describes the score. The score units are
completely unspecified, but for sequence similarities, it is typically percent
identity. Annotations that do not have a score can use "."

strand

For those annotations which are strand-specific, this field is the strand on
which the annotation resides. It is "+" for the forward strand, "-" for the
reverse strand, or "." for annotations that are not stranded.

phase

For annotations that are linked to proteins, this field describes the phase of
the annotation on the codons. It is a number from 0 to 2, or "." for features
that have no phase.

group

GFF provides a simple way of generating annotation hierarchies ("is composed
of" relationships) by providing a group field. The group field contains the
class and ID of an annotation which is the logical parent of the current one.
In the example given above, the group is the Transcript named "R119.7".
The group field is also used to store information about the target of sequence
similarity hits, and miscellaneous notes. See the next section for a
description of how to describe similarity targets.

The sequences used to establish the coordinate system for annotations can
correspond to sequenced clones, clone fragments, contigs or super-contigs.

In addition to a group ID, the GFF format allows annotations to have a group
class. This makes sure that all groups are unique even if they happen to share
the same name. For example, you can have a GenBank accession named AP001234 and
a clone named AP001234 and distinguish between them by giving the first one a
class of Accession and the second a class of Clone.

You should use double-quotes around the group name or class if it contains
white space.


https://uswest.ensembl.org/info/website/upload/gff.html
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GFF (General Feature Format) format consists of one line per feature, each
containing 9 columns of data, plus optional track definition lines. The
following documentation is based on the Version 2 specifications.

The GTF (General Transfer Format) is identical to GFF version 2.

Fields
Track lines
More information
Fields

Fields must be tab-separated. Also, all but the final field in each feature
line must contain a value; "empty" columns should be denoted with a '.'

seqname - name of the chromosome or scaffold; chromosome names can be given
with or without the 'chr' prefix. Important note: the seqname must be one used
within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such
as a scaffold ID, without any additional content such as species or assembly.
See the example GFF output below.

source - name of the program that generated this feature, or the data source
(database or project name)

feature - feature type name, e.g. Gene, Variation, Similarity

start - Start position of the feature, with sequence numbering starting at 1.

end - End position of the feature, with sequence numbering starting at 1.

score - A floating point value.

strand - defined as + (forward) or - (reverse).

frame - One of '0', '1' or '2'. '0' indicates that the first base of the
feature is the first base of a codon, '1' that the second base is the first
base of a codon, and so on..

attribute - A semicolon-separated list of tag-value pairs, providing additional
information about each feature.


https://genome.ucsc.edu/FAQ/FAQformat.html#format3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GFF (General Feature Format) lines are based on the Sanger GFF2 specification.
GFF lines have nine required fields that must be tab-separated. If the fields
are separated by spaces instead of tabs, the track will not display correctly.
For more information on GFF format, refer to Sanger's GFF page.

Note that there is also a GFF3 specification that is not currently supported by
the Browser. All GFF tracks must be formatted according to Sanger's GFF2
specification.
If you would like to obtain browser data in GFF (GTF) format, please refer to
Genes in gtf or gff format on the Wiki.

Here is a brief description of the GFF fields:

seqname - The name of the sequence. Must be a chromosome or scaffold.

source - The program that generated this feature.

feature - The name of this type of feature. Some examples of standard feature
types are "CDS" "start_codon" "stop_codon" and "exon"li>

start - The starting position of the feature in the sequence. The first base is
numbered 1.

end - The ending position of the feature (inclusive).

score - A score between 0 and 1000. If the track line useScore attribute is set
to 1 for this annotation data set, the score value will determine the level of
gray in which this feature is displayed (higher numbers = darker gray). If
there is no score value, enter ":.":.

strand - Valid entries include "+", "-", or "." (for don't know/don't care).

frame - If the feature is a coding exon, frame should be a number between 0-2
that represents the reading frame of the first base. If the feature is not a
coding exon, the value should be ".".

group - All lines with the same group are linked together into a single item.

"""

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

