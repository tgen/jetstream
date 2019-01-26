import os
import tempfile
import jetstream
from jetstream.cli.jetstream import main

from unittest import TestCase

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestCliExt(TestCase):
    """Tests that run workflow templates stored externally in
    test/test_templates """

    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(TestCliExt, self).setUp()
        self.templates_dir = os.path.abspath(os.path.join(TESTS_DIR, 'test_templates'))
        self.variables = os.path.abspath(os.path.join(TESTS_DIR, 'test_templates', 'variables.yaml'))
        self._original_dir = os.getcwd()
        self._temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self._temp_dir.name)

    def tearDown(self):
        os.chdir(self._original_dir)
        self._temp_dir.cleanup()

    def test_should_pass(self):
        templates_dir = os.path.join(self.templates_dir, 'should_pass')
        templates = os.listdir(templates_dir)

        for f in templates:
            t = os.path.join(templates_dir, f)

            with self.subTest(msg=t):
                main(['run', t, '--variables', self.variables, '--backend', 'local'])



class TestCli(TestCase):
    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(TestCli, self).setUp()
        self._original_dir = os.getcwd()
        self._temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self._temp_dir.name)

    def tearDown(self):
        os.chdir(self._original_dir)
        self._temp_dir.cleanup()

    def test_init(self):
        init.main()

    def test_run(self):
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: hostname\n')

        run.main([
            'testwf.jst',
            '--backend', 'local'
        ])

    def test_run_w_vars(self):
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: echo {{ name }}\n  stdout: /dev/null')

        run.main([
            'testwf.jst',
            '--backend', 'local',
            '--str:name', 'Philip J. Fry',
            '--bool:ok', 'true',
            '--int:number', '42',
            '--float:number2', '3.14'
        ])

    def test_project_tasks(self):
        jetstream.Project(new=True)

        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: hostname\n')

        run.main([
            'testwf.jst',
            '--backend', 'local'
        ])

        project.main(['tasks'])
