""" This should bridge the workflow components and the execution environment. 
Module names should be linked to the command that executes them. We need to decide
on exactly how the modules behave:

    - Do they receive a path to the project as an argument?
    - Can they assume that they will always start in a project directory?
    - Should they read project data from some sort of manifest?
    - Or, is piped to stdin?
    
"""    
import jetstream.settings as settings
from random import random


def faux_cmd(msg):
    delay = int(random() * 45)
    cmd = ['bash', '-c', 'sleep {} && echo hello from {} $$'.format(delay, msg)]
    return cmd


def lookup_module_cmd(name):
    modules = {
        'dna_align': {
            'cmd': ['bash', '-c', 'sleep 30 && echo dna_align done']
        },
        'joint_indel_realign': {
            'cmd': ['bash', '-c', 'sleep 30 && echo joint_indel_realign done']
        },
        'germline_variant_calling': {
            'cmd': ['bash', '-c', 'sleep 30 && echo germline_variant_calling done']
        },
        'somatic_variant_calling': {
            'cmd': ['bash', '-c', 'sleep 30 && echo somatic_variant_calling done']
        },
        'snpeff': {
            'cmd': ['bash', '-c', 'sleep 30 && echo snpeff done']
        },
        'rna_quant_htseq': {
            'cmd': ['bash', '-c', 'sleep 30 && echo rna_quant_htseq done']
        },
        'multiqc': {
            'cmd': ['bash', '-c', 'sleep 30 && echo multiqc done']
        },
        'rna_align_star': {
            'cmd': ['bash', '-c', 'sleep 30 && echo rna_align_star done']
        }
    }

    return modules[name]

