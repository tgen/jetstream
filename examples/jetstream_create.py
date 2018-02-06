#!/usr/bin/env python3
import logging
from jetstream.workflow import Workflow, serialize_json

# TODO These should be written with a unittest framework

def test():
    wf = Workflow()
    bwa_node = wf.add_component('dna_alignment/bwa_mem.yaml')
    md_node = wf.add_component_after('dna_alignment/mark_duplicates.yaml', bwa_node)
    wf.add_component_after('dna_alignment/bqsr.yaml', md_node)
    wf.add_component('rna_alignment/star_alignment.yaml')
    return wf


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    log.critical('Logging started')
    wf = test()
    print(serialize_json(wf))

