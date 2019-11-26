import os
import tempfile
import jetstream
from unittest import TestCase
from jetstream.cli import main as cli_main

jetstream.settings.clear()
jetstream.settings.read(user=False)
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_TEMPLATES = os.path.join(TESTS_DIR, 'templates')


class TestCliRunTemplates(TestCase):
    """Tests that run workflow templates stored externally in
    test/test_templates """
    longMessage = True

    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(TestCliRunTemplates, self).setUp()
        self.original_dir = os.getcwd()
        self.temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.temp_dir.name)

    def tearDown(self):
        os.chdir(self.original_dir)
        self.temp_dir.cleanup()

    def template(self, template_filename):
        """Runs a template from the test template directory"""
        path = os.path.join(TEST_TEMPLATES, template_filename)
        args = ('run', path)
        cli_main(args)
       
    def test_helloworld_1(self):
        """single helloworld task"""
        self.template('helloworld_1.jst')

    def test_dependencies_1(self):
        """dependencies can be declared with before/after"""
        self.template('dependencies_1.jst')
    
    def test_dependencies_2(self):
        """dependencies can be declared with input/output"""
        self.template('dependencies_2.jst')

    def test_dependencies_3(self):
        """tasks with dependencies can be added during run with exec"""
        self.template('dependencies_3.jst')

    def test_inheritance_1(self):
        """templates can include code from other templates"""
        self.template('inheritance_1.jst')

    def test_logging_1(self):
        """templates can log to stderr with log global fn"""
        self.template('logging_1.jst')

    def test_mapping_1(self):
        """templates can include properties mapping at the top"""
        self.template('mapping_1.jst')

    def test_retry_1(self):
        """tasks can include retry directive that allows tasks to fail and
        then be run again"""
        self.template('retry_1.jst')

    def test_stress_1(self):
        """runs should not crash due to forking limits"""
        self.template('stress_1.jst')

    def test_stress_2(self):
        """runs should be able to process lots of concurrent tasks"""
        self.template('stress_2.jst')

    def test_stress_3(self):
        """tasks with no command should complete very fast"""
        self.template('stress_3.jst')
