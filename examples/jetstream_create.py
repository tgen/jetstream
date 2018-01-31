#!/usr/bin/env python3
import logging
from jetstream.workflow import Workflow

# TODO These should be written with a unittest framework

def test():
    wf = Workflow('test_wf')
    wf.add_module('pegasusPipe/pegasus_dnaAlign.sh')
    wf.add_module_after('pegasusPipe/pegasus_dnaAlign.sh', 'pegasusPipe/pegasus_indelRealign.sh')
    wf.add_module_after('pegasusPipe/pegasus_dnaAlign.sh', 'pegasusPipe/pegasus_haplotypeCaller.sh')
    wf.add_module_after('pegasusPipe/pegasus_indelRealign.sh', 'pegasusPipe/pegasus_mutect.sh')
    wf.add_module_after('pegasusPipe/pegasus_mutect.sh', 'pegasusPipe/pegasus_vcfMerger.sh')
    wf.add_module_before('pegasusPipe/pegasus_vcfMerger.sh', 'pegasusPipe/pegasus_haplotypeCaller.sh')
    wf.add_module('pegasusPipe/pegasus_htSeq.sh')
    wf.add_module_before('pegasusPipe/pegasus_htSeq.sh', 'pegasusPipe/pegasus_rnaAlign.sh')
    wf.add_module('pegasusPipe/pegasus_summaryStats.sh')
    wf.add_module_before(
        'pegasusPipe/pegasus_summaryStats.sh',
        'pegasusPipe/pegasus_vcfMerger.sh',
        'pegasusPipe/pegasus_htSeq.sh'
    )
    return wf


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    log.critical('Logging started')
    wf = test()
    print(wf.serialize_pydot())

