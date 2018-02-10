#!/usr/bin/env python3
import unittest
import jetstream


class TestWorkflowModule(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestWorkflowModule, cls).setUpClass()
        try:
            jetstream.plugins.clone('https://github.com/ryanrichholt/test_components.git')
        except Exception:
            jetstream.plugins.update()

    def test_create_workflow(self):
        wf = jetstream.Workflow()

    def test_add_component(self):
        wf = jetstream.Workflow()
        bwa_node = wf.add_component('test_components/bwa_mem.yaml')

    def test_add_component_after(self):
        wf = jetstream.Workflow()
        bwa_node = wf.add_component('test_components/bwa_mem.yaml')
        md_node = wf.add_component_after('test_components/mark_duplicates.yaml', bwa_node)

    def test_add_component_before(self):
        wf = jetstream.Workflow()
        bwa_node = wf.add_component('test_components/bwa_mem.yaml')
        md_node = wf.add_component_before('test_components/mark_duplicates.yaml', bwa_node)
