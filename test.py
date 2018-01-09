import logging
from jetstream.workflow import Workflow
from jetstream.async_runner import run


def test():
    wf = Workflow('test_wf')
    wf.add_module('dna_align')
    wf.add_module_after('dna_align', 'joint_indel_realign')
    wf.add_module_after('dna_align', 'germline_variant_calling')
    wf.add_module_after('joint_indel_realign', 'somatic_variant_calling')
    wf.add_module_after('somatic_variant_calling', 'snpeff')
    wf.add_module_before('snpeff', 'germline_variant_calling')
    wf.add_module_before('rna_quant_htseq', 'rna_align_star')
    wf.add_module('multiqc')
    wf.add_module_before('multiqc', 'snpeff', 'rna_quant_htseq')
    return wf


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    log.critical('Logging started')
    wf = test()
    print(wf)
    # run(wf, debug=False)
