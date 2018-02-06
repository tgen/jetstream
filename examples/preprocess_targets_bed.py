#!/usr/bin/env python3
import sys
import argparse
from jetstream.formats import intervals
from jetstream.utils import remove_prefix


grcm38_mm10_UCSC_TO_ENSEMBL = {
    'chr1': '1',
    'chr2': '2',
    'chr3': '3',
    'chr4': '4',
    'chr5': '5',
    'chr6': '6',
    'chr7': '7',
    'chr8': '8',
    'chr9': '9',
    'chr10': '10',
    'chr11': '11',
    'chr12': '12',
    'chr13': '13',
    'chr14': '14',
    'chr15': '15',
    'chr16': '16',
    'chr17': '17',
    'chr18': '18',
    'chr19': '19',
    'chrX': 'X',
    'chrY': 'Y',
    'chrM': 'MT',
    'chr1_GL456210_random': 'GL456210.1',
    'chr1_GL456211_random': 'GL456211.1',
    'chr1_Gl456212_random': 'GL456212.1',
    'chr1_GL456213_random': 'GL456213.1',
    'chr1_GL456221_random': 'GL456221.1',
    'chr4_GL456216_random': 'GL456216.1',
    'chr4_GL456350_random': 'GL456350.1',
    'chr4_JH584292_random': 'JH584292.1',
    'chr4_JH584293_random': 'JH584293.1',
    'chr4_JH584294_random': 'JH584294.1',
    'chr4_JH584295_random': 'JH584295.1',
    'chr5_GL456354_random': 'GL456354.1',
    'chr5_JH584296_random': 'GH584296.1',
    'chr5_JH584297_random': 'GH584297.1',
    'chr5_JH584298_random': 'GH584298.1',
    'chr5_JH584299_random': 'GH584299.1',
    'chr7_GL456219_random': 'GL456219.1',
    'chrUn_GL456239': 'GL456239.1',
    'chrUn_GL456359': 'GL456359.1',
    'chrUn_GL456360': 'GL456360.1',
    'chrUn_GL456366': 'GL456366.1',
    'chrUn_GL456367': 'GL456367.1',
    'chrUn_GL456368': 'GL456368.1',
    'chrUn_GL456370': 'GL456370.1',
    'chrUn_GL456372': 'GL456372.1',
    'chrUn_GL456378': 'GL456378.1',
    'chrUn_GL456379': 'GL456379.1',
    'chrUn_GL456381': 'GL456381.1',
    'chrUn_GL456382': 'GL456382.1',
    'chrUn_GL456383': 'GL456383.1',
    'chrUn_GL456385': 'GL456385.1',
    'chrUn_GL456387': 'GL456387.1',
    'chrUn_GL456389': 'GL456389.1',
    'chrUn_GL456390': 'GL456390.1',
    'chrUn_GL456392': 'GL456392.1',
    'chrUn_GL456393': 'GL456393.1',
    'chrUn_GL456394': 'GL456394.1',
    'chrUn_GL456396': 'GL456396.1',
    'chrUn_JH584304': 'JH584304.1',
    'chrX_GL456233_random': 'GL456233.1',
    'chrY_JH584300_random': 'JH584300.1',
    'chrY_JH584301_random': 'JH584301.1',
    'chrY_JH584302_random': 'JH584302.1',
    'chrY_JH584303_random': 'JH584303.1',
}

grch37_hg19_UCSC_TO_ENSEMBL = {
    'chr1': '1',
    'chr2': '2',
    'chr3': '3',
    'chr4': '4',
    'chr5': '5',
    'chr6': '6',
    'chr7': '7',
    'chr8': '8',
    'chr9': '9',
    'chr10': '10',
    'chr11': '11',
    'chr12': '12',
    'chr13': '13',
    'chr14': '14',
    'chr15': '15',
    'chr16': '16',
    'chr17': '17',
    'chr18': '18',
    'chr19': '19',
    'chr20': '20',
    'chr21': '21',
    'chr22': '22',
    'chrX': 'X',
    'chrY': 'Y',
    'chrM': 'MT',
    'chr1_gl000191_random': 'GL000191.1',
    'chr1_gl000192_random': 'GL000192.1',
    'chr4_gl000193_random': 'GL000193.1',
    'chr4_gl000194_random': 'GL000194.1',
    'chr7_gl000195_random': 'GL000195.1',
    'chr8_gl000196_random': 'GL000196.1',
    'chr8_gl000197_random': 'GL000197.1',
    'chr9_gl000198_random': 'GL000198.1',
    'chr9_gl000199_random': 'GL000199.1',
    'chr9_gl000200_random': 'GL000200.1',
    'chr9_gl000201_random': 'GL000201.1',
    'chr11_gl000202_random': 'GL000202.1',
    'chr17_gl000203_random': 'GL000203.1',
    'chr17_gl000204_random': 'GL000204.1',
    'chr17_gl000205_random': 'GL000205.1',
    'chr17_gl000206_random': 'GL000206.1',
    'chr18_gl000207_random': 'GL000207.1',
    'chr19_gl000208_random': 'GL000208.1',
    'chr19_gl000209_random': 'GL000209.1',
    'chr21_gl000210_random': 'GL000210.1',
    'chrUn_gl000211': 'GL000211.1',
    'chrUn_gl000212': 'GL000212.1',
    'chrUn_gl000213': 'GL000213.1',
    'chrUn_gl000214': 'GL000214.1',
    'chrUn_gl000215': 'GL000215.1',
    'chrUn_gl000216': 'GL000216.1',
    'chrUn_gl000217': 'GL000217.1',
    'chrUn_gl000218': 'GL000218.1',
    'chrUn_gl000219': 'GL000219.1',
    'chrUn_gl000220': 'GL000220.1',
    'chrUn_gl000221': 'GL000221.1',
    'chrUn_gl000222': 'GL000222.1',
    'chrUn_gl000223': 'GL000223.1',
    'chrUn_gl000224': 'GL000224.1',
    'chrUn_gl000225': 'GL000225.1',
    'chrUn_gl000226': 'GL000226.1',
    'chrUn_gl000227': 'GL000227.1',
    'chrUn_gl000228': 'GL000228.1',
    'chrUn_gl000229': 'GL000229.1',
    'chrUn_gl000230': 'GL000230.1',
    'chrUn_gl000231': 'GL000231.1',
    'chrUn_gl000232': 'GL000232.1',
    'chrUn_gl000233': 'GL000233.1',
    'chrUn_gl000234': 'GL000234.1',
    'chrUn_gl000235': 'GL000235.1',
    'chrUn_gl000236': 'GL000236.1',
    'chrUn_gl000237': 'GL000237.1',
    'chrUn_gl000238': 'GL000238.1',
    'chrUn_gl000239': 'GL000239.1',
    'chrUn_gl000240': 'GL000240.1',
    'chrUn_gl000241': 'GL000241.1',
    'chrUn_gl000242': 'GL000242.1',
    'chrUn_gl000243': 'GL000243.1',
    'chrUn_gl000244': 'GL000244.1',
    'chrUn_gl000245': 'GL000245.1',
    'chrUn_gl000246': 'GL000246.1',
    'chrUn_gl000247': 'GL000247.1',
    'chrUn_gl000248': 'GL000248.1',
    'chrUn_gl000249': 'GL000249.1',
    'chr4_ctg9_hap1': 'HSCHR4_1',
    'chr6_apd_hap1': 'HSCHR6_MHC_APD',
    'chr6_cox_hap2': 'HSCHR6_MHC_COX',
    'chr6_dbb_hap3': 'HSCHR6_MHC_DBB',
    'chr6_mann_hap4': 'HSCHR6_MHC_MANN',
    'chr6_mcf_hap5': 'HSCHR6_MHC_MCF',
    'chr6_qbl_hap6': 'HSCHR6_MHC_QBL',
    'chr6_ssto_hap7': 'HSCHR6_MHC_SSTO',
    'chr17_ctg5_hap1': 'HSCHR17_1',
}

def preprocess_bed(path, replace_underscores=False, replace_mt=False,
                   ucsc_to_ensembl=True):
    """ Preprocessing used for canine bed files from Agilent """
    bed = intervals.read_bed(path)

    # Iterate over every interval in the bed file
    for i in bed:
        if ucsc_to_ensembl:
            i['seqname'] = grch37_hg19_UCSC_TO_ENSEMBL.get(i['seqname'], i['seqname'])
        if replace_underscores:
            i['seqname'] = i['seqname'].replace('_', '.')
        if replace_mt:
            if i['seqname'] == 'M':
                i['seqname'] = 'MT'

    print(intervals.to_bed(bed))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Preprocess bed files prior to exome kit refpack creation. "
                    "This will remove 'chr' prefixes from chromosome names and "
                    "replace chromosome 'M' with 'MT'"
    )

    parser.add_argument('bed', help='path to a bed file to preprocess')

    parser.add_argument('--replace-underscores', action='store_true', default=False,
                        help="Replace underscores with dots. Useful for some "
                             "Agilent files that don't match standard chrom "
                             "names, but can cause issues on most others.")

    parser.add_argument('--replace-mt', action='store_true', default=False,
                        help="Replace chromosome name 'M' with 'MT' to match "
                             "reference genomes.")

    parser.add_argument('--ucsc-to-ensembl',
                        action='store_true', default=False,
                        help="replace UCSC sequence names with Ensembl.")

    args = parser.parse_args()

    preprocess_bed(
        path=args.bed,
        replace_underscores=args.replace_underscores,
        replace_mt=args.replace_mt,
        ucsc_to_ensembl=args.ucsc_to_ensembl,
    )



