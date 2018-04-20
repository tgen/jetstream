import json
import logging
import re
import sys

log = logging.getLogger(__name__)

CMD_TEMPLATE = r"""# Align with BWA mem and sort with samtools

bwa mem -R "{READGROUP}" -M -t8 {REF} {R1} {R2} | \
  samtools sort -O BAM -o {OUTPATH} -
  
"""

def bwa(r1fastq, r2fastq, reference_index, read_group, outpath):
    final_cmd = CMD_TEMPLATE.format(
        READGROUP=read_group,
        REF=reference_index,
        R1=r1fastq,
        R2=r2fastq,
        OUTPATH=outpath
    )

    j = utils.batch_schedulers.slurm.easy(
        final_cmd,
        '-c8', '--mem3000',
        module_load='bwa/0.7.12 samtools/1.4.1'
    )
    return j.wait()


def launch_bwa_on_fastq_record(fastq, reference_index):
    """Resolves necessary input arguments from a fastq record
    and launches a bwa job"""

    read_group = utils.read_group(
        ID=fastq['RG_ID'],
        CN=fastq['RG_CN'],
        KS=fastq['RG_KS'],
        SM=fastq['RG_SM'],
        LB=fastq['RG_LB'],
        PU=fastq['RG_PU']
    )

    r1fastq = fastq['path']
    r2fastq = re.subn('R1_001', 'R2_001', r1fastq)

    outpath="{sm}/{sm}.{lb}.bam".format(sm=fastq['RG_SM'], lb=fastq['RG_LB'])

    return bwa(r1fastq, r2fastq, reference_index, read_group, outpath)


def main():
    json_data = sys.argv[1]
    reference_index = sys.argv[2]

    print(json.loads(json_data))


if __name__ == '__main__':
    main()
