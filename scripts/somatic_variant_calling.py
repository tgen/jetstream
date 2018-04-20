import argparse

from jetstream import Project
from jetstream.utils import slurm

MUTECT2_CMD = r"""
set -eu

gatk Mutect2 \
   -nct 8 \
   -R {reference_fa} \
   -I {tumor_bam} \
   -tumor {tumor_sample_name} \
   -I {normal_bam} \
   -normal {normal_sample_name} \
   --germline-resource {af_only_gnomad_vcf_gz} \
   --af-of-alleles-not-in-resource 0.00003125 \
   --panel-of-normals {pon_vcf_gz} \
   --intervals {interval_list} \
   --verbosity INFO \
   -O {out_path}

"""


def mutect2(tumor_bam, normal_bam):
    final_cmd =MUTECT2_CMD.format(
        #TODO
    )

    j = slurm.easy(
        final_cmd,
        '-c8',
        module_load='gatk/4.0.1.2'
    )

    return j


def launch_on_all_pairs_in_project():
    slurm_jobs = []
    project = Project()

    reference = project.reference('reference_genome_bwa_index')
    known_sites = project.reference('known_sites')

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
        description="Launch Somatic Variant Calling on a Jetstream project"
    )

    parser.add_argument('--overwrite', action='store_true', default=False)

    return parser


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)


if __name__ == '__main__':
    main()
