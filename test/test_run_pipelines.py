import os
import tempfile
from unittest import TestCase
import jetstream
from jetstream.cli import main as cli_main

jetstream.settings.clear()
jetstream.settings.read(user=False)
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PIPELINES = os.path.join(TESTS_DIR, 'pipelines')


class TestCliRunPipelines(TestCase):
    """Tests that run workflow templates stored externally in
    test/test_templates """
    longMessage = True

    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(TestCliRunPipelines, self).setUp()
        self.original_dir = os.getcwd()
        self.temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.temp_dir.name)

    def tearDown(self):
        os.chdir(self.original_dir)
        self.temp_dir.cleanup()

    def pipeline(self, pipeline_name):
        jetstream.settings['pipelines']['searchpath'] = TEST_PIPELINES
        args = ('pipelines', pipeline_name)
        cli_main(args)

    def test_foopipe_1(self):
        """the most basic pipeline requires pipeline.yaml and main template"""
        self.pipeline('foopipe_1')

    def test_foopipe_2(self):
        """other pipeline features include bin directory"""
        self.pipeline('foopipe_2')

    def test_foopipe_3(self):
        """pipelines nested in other pipelines can be discovered as well"""
        # barpipe_1 is nested in foopipe_3
        self.pipeline('barpipe_1')

    def test_foopipe_4(self):
        """pipeline variables are exported and allow files to be accessed"""
        self.pipeline('foopipe_4')
