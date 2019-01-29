import os
import tempfile
import jetstream
from unittest import TestCase
from jetstream.cli.jetstream import main as cli_main

CMD_ARGS = ['--logging', 'silent']
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_TEMPLATES = os.path.join(TESTS_DIR, 'templates')
TEST_VARIABLES = os.path.join(TESTS_DIR, 'templates', 'variables.yaml')


class TestCliExt(TestCase):
    """Tests that run workflow templates stored externally in
    test/test_templates """

    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(TestCliExt, self).setUp()
        self.original_dir = os.getcwd()
        self.temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.temp_dir.name)

    def tearDown(self):
        os.chdir(self.original_dir)
        self.temp_dir.cleanup()

    def test_should_pass(self):
        templates_dir = os.path.join(TEST_TEMPLATES, 'should_pass')
        templates = os.listdir(templates_dir)

        for f in templates:
            t = os.path.join(templates_dir, f)

            with self.subTest(msg=t):
                args = CMD_ARGS + [
                    'run',
                    t,
                    '--variables', TEST_VARIABLES
                ]

                cli_main(args)


class TestCli(TestCase):
    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(TestCli, self).setUp()
        self.original_dir = os.getcwd()
        self.temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.temp_dir.name)

    def tearDown(self):
        os.chdir(self.original_dir)
        self.temp_dir.cleanup()

    def test_init(self):
        args = CMD_ARGS +[
            'init',
        ]

        cli_main(args)

    def test_run(self):
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: hostname\n')

        args = CMD_ARGS + [
            'run',
            'testwf.jst',
        ]

        cli_main(args)

    def test_run_w_vars(self):
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: echo {{ name }}\n  stdout: /dev/null')

        args = CMD_ARGS + [
            'run',
            'testwf.jst',
            '--str:name', 'Philip J. Fry',
            '--bool:ok', 'true',
            '--int:number', '42',
            '--float:number2', '3.14'
        ]

        cli_main(args)

    def test_project_tasks(self):
        jetstream.new_project()

        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: hostname\n')

        args = CMD_ARGS + [
            'run',
            'testwf.jst',
        ]

        cli_main(args)

        args = CMD_ARGS + [
            'project',
            'tasks',
        ]

        cli_main(args)
