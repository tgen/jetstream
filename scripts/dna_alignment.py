#!/usr/bin/env python3
"""Generate preliminary BAMs from FASTQs. When executed from the
command line, this module will launch Slurm jobs for every fastq
set in the current project. It also contains functions for launching
individual jobs that can be used via import.

Requires a reference data set in the project with the following:

TODO

"""
import argparse
import logging
import os
import re

from jetstream import Project
from jetstream import read_group
from jetstream.utils.batch_schedulers import slurm

log = logging.getLogger(__name__)

ALIGNMENT_CMD = r"""

bwa mem -R "{RGTAG}" -M -t8 {REF} {R1} {R2} | samtools sort -O BAM -o {OUT_PATH} - 

gatk BaseRecalibrator \
    --reference {REF} \
    --known-sites {KNOWN_SITES} \
    --input {OUT_PATH} \
    --output {OUT_PATH}.recal_data.table

gatk ApplyBQSR \
   --reference {REF} \
   --input {OUT_PATH} \
   --bqsr-recal-file {OUT_PATH}.recal_data.table \
   -O {OUT_PATH}.bam

"""

def align_and_recalibrate(r1_fastq, r2_fastq, rg, out, reference, known_sites):
    """Launch bwa mem on a set of fastqs using slurm.
    :param r1_fastq: path to read 1 fastq
    :param r2_fastq: path to read 2 fastq
    :param rg: read group dictionary
    :param out: out path
    :param reference: path to bwa reference index
    :return: SlurmJob
    """
    # Builds a read group line from a dictionary containing known tags
    rg_tag = read_group(**rg)

    final_cmd = ALIGNMENT_CMD.format(
        RGTAG=rg_tag,
        REF=reference,
        R1=r1_fastq,
        R2=r2_fastq,
        KNOWN_SITES=known_sites,
        OUT_PATH=out
    )

    j = slurm.easy(
        final_cmd,
        '-c8',
        module_load='bwa/0.7.12 samtools/1.4.1 gatk/4.0.1.2'
    )

    return j


def launch_bwamem_on_data(data, reference, known_sites, overwrite):
    r1_fastq = data['path']
    r2_fastq, n = re.subn('R1_001', 'R2_001', r1_fastq)

    if n != 1:
        raise ValueError(
            'Fastq read name substitution failed for {}'.format(r1_fastq))

    os.makedirs('{sample_name}/{assay}'.format(
        sample_name=data['sample_name'],
        assay=data['assay']
    ), exist_ok=True)

    out_path = '{sample_name}/{assay}/{sample_name}.bam'.format(
        sample_name=data['sample_name'],
        assay=data['assay']
    )

    if os.path.exists(out_path) and not overwrite:
        raise FileExistsError(out_path)

    job = align_and_recalibrate(
        r1_fastq=r1_fastq,
        r2_fastq=r2_fastq,
        rg=data['read_group'],
        reference=reference,
        known_sites=known_sites,
        out=out_path
    )

    return job


def launch_alignment_on_project(overwrite):
    """Launches bwa mem on every data applicable item in a project"""
    slurm_jobs = []
    project = Project()
    reference = project['reference']['reference_genome_bwa_index']
    known_sites = project['reference']['known_sites']

    to_align = project.data(type='FQ', read_style='paired-end')
    if not to_align:
        raise RuntimeError('No data to align found in project')

    for data in to_align:
        try:
            job = launch_bwamem_on_data(data, reference, known_sites, overwrite)
            slurm_jobs.append(job)
        except Exception as e:
            log.exception(e)

    return slurm_jobs


def arg_parser():
    parser = argparse.ArgumentParser(
        description="This launches DNA Alignnment on a Jetstream project"
    )

    parser.add_argument('--overwrite', action='store_true', default=False)

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)

    jobs = launch_alignment_on_project(overwrite=args.overwrite)
    print('Waiting on slurm jobs {}'.format(jobs))
    slurm.wait_for(*jobs)
    print('Done with bwa mem!')

    # launch merge on each sample

    # launch indel realign on each sample

    # launch mark duplicates on each sample


if __name__ == '__main__':
    main()
