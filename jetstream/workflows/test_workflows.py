from jetstream import workflows


def somatic_variant_calling():
    wf = workflows.Workflow(name="Somatic Variant Calling")
    wf.add_node("Mark Duplicates", 'test_components.git/mark_duplicates.yaml')
    wf.add_node_after("Mutect", 'test_components.git/mutect.yaml',
                      "Mark Duplicates")
    wf.add_node_after("Strelka", 'test_components.git/strelka.yaml',
                      "Mark Duplicates")
    wf.add_node_after("Seurat", 'test_components.git/seurat.yaml',
                      "Mark Duplicates")
    wf.add_node_after("Merge VCFs", 'test_components.git/vcf_merger.yaml',
                      "Mutect", "Strelka", "Seurat")
    return wf


def dna_alignments():
    wf = workflows.Workflow(name="DNA Alignment")
    wf.add_node("BWA-MEM", 'test_components.git/bwa_mem.yaml')
    wf.add_node_after("GATK BaseRecalibrator",
                      'test_components.git/recalibrate.yaml',
                      "BWA-MEM")
    wf.add_node_after("Picard Mark Duplicates",
                      'test_components.git/mark_duplicates.yaml',
                      "GATK BaseRecalibrator")
    return wf
