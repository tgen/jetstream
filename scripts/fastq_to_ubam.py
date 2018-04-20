#!/usr/bin/env python3
"""Converts all paired-end fastqs in a project to unmapped bams"""
import jetstream


def fastq_to_sam(data):
    cmd_template = r"""
    gatk FastqToSam \
        -F1 {r1fq} \
        -F2 {r2fq} \
        -O  {sample_name}/{sample_name}.bam \
        -RG {rg_id} \
        -SM {rg_sn} \
        -LB {rg_lb} \
        -PU {rg_pu} \
        -PL {rg_pl} \
        -CN {rg_cn}
    """
    cmd = cmd_template.format(
        r1fq=data['path'],
        r2fq=data['path'].replace('R1_001.fastq', 'R2_001.fastq'),
        sample_name=data['sample_name'],
        rg_id=data['read_group']['ID'],
        rg_sn=data['read_group']['SM'],
        rg_lb=data['read_group']['LB'],
        rg_pu=data['read_group']['PU'],
        rg_pl=data['read_group']['PL'],
        rg_cn=data['read_group']['CN']
    )

    j = jetstream.utils.slurm.easy(
        cmd, '--mem', '24000', module_load='gatk/4.0.1.2')
    print(j)


if __name__ == "__main__":
    p = jetstream.Project()
    for data in p.list_data(type='FQ', read_style='paired_end'):
        fastq_to_sam(data)
