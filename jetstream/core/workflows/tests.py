#!/usr/bin/env python3
import os
import tempfile
import unittest

from jetstream import plugins
from jetstream.core import Project, init, workflows


class TestWorkflowModule(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestWorkflowModule, cls).setUpClass()
        try:
            plugins.clone('https://github.com/ryanrichholt/test_components.git')
        except Exception:
            plugins.update()

    def test_create_workflow(self):
        workflows.Workflow()

    def test_add_component(self):
        wf = workflows.Workflow()
        wf.add_node('test_components.git/bwa_mem.yaml')

    def test_add_component_after(self):
        wf = workflows.Workflow()
        bwa_node = wf.add_node('test_components.git/bwa_mem.yaml')
        wf.add_node_after('test_components.git/mark_duplicates.yaml', bwa_node)

    def test_add_component_before(self):
        wf = workflows.Workflow()
        bwa_node = wf.add_node('test_components.git/bwa_mem.yaml')
        wf.add_node_before('test_components.git/mark_duplicates.yaml', bwa_node)


class TestWorkflowRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestWorkflowRun, cls).setUpClass()

        # Make sure we have the latest test components
        try:
            plugins.clone('https://github.com/ryanrichholt/test_components.git')
        except Exception:
            plugins.update()

    def test_wf_run(self):
        tmp_dir = tempfile.TemporaryDirectory()
        print(tmp_dir.name)
        os.chdir(tmp_dir.name)

        workflow = workflows.Workflow()
        node1 = workflow.add_node('test_components.git/bwa_mem.yaml')
        node2 = workflow.add_node_after('test_components.git/bwa_mem.yaml', node1)
        workflow.add_node_after('test_components.git/bwa_mem.yaml', node2)

        init()
        p = Project()
        p.run(workflow)
        p.latest_run()
        tmp_dir.cleanup()

