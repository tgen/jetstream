import os
import sys
import tempfile
import jetstream
from jetstream.cli.subcommands import init, run, pipelines, project

from unittest import TestCase


class TestCliExt(TestCase):
    """Tests that run workflow templates stored externally in
    test/test_templates """
    templates_dir = os.path.join('test', 'test_templates')
    templates_config = os.path.join('test', 'test_templates', 'config.yaml')

    def test_should_pass(self):
        templates_dir = os.path.join(self.templates_dir, 'should_pass')
        templates = os.listdir(templates_dir)

        for f in templates:
            template = os.path.join(templates_dir, f)

            with self.subTest(msg=template):
                with self.assertRaises(SystemExit) as cm:
                    run.main([
                        template,
                        '--config', self.templates_config
                    ])

                self.assertEqual(cm.exception.code, 0)

    def test_should_fail(self):
        templates_dir = os.path.join(self.templates_dir, 'should_fail')
        templates = os.listdir(templates_dir)

        for f in templates:
            template = os.path.join(templates_dir, f)

            with self.subTest(msg=template):
                with self.assertRaises(Exception) as cm:
                    run.main([
                        template,
                        '--config', self.templates_config
                    ])



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

        with self.assertRaises(SystemExit) as cm:
            run.main([
                'testwf.jst',
                '--backend',
                'local'
            ])

        self.assertEqual(cm.exception.code, 0)

    def test_run_w_vars(self):
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: echo {{ name }}\n  stdout: /dev/null')

        with self.assertRaises(SystemExit) as cm:
            run.main([
                'testwf.jst',
                '--backend',
                'local',
                '--str:name',
                'Philip J. Fry'
            ])

        self.assertEqual(cm.exception.code, 0)

    def test_pipelines(self):
        jetstream.Project(new=True)

        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: hostname\n')

        with self.assertRaises(SystemExit) as cm:
            pipelines.main([
                'testwf.jst',
                '--backend',
                'local'
            ])

        self.assertEqual(cm.exception.code, 0)

    def test_project_tasks(self):
        p = jetstream.Project(new=True)

        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: hostname\n')

        with self.assertRaises(SystemExit) as cm:
            pipelines.main([
                'testwf.jst',
                '--backend',
                'local'
            ])

        self.assertEqual(cm.exception.code, 0)

        project.main(['tasks'])
