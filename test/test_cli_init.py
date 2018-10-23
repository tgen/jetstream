import os
import sys
import tempfile
import jetstream
from jetstream.cli.subcommands import init, run, pipelines, project

from unittest import TestCase


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

        with self.assertRaises(SystemExit):
            run.main([
                'testwf.jst',
                '--backend',
                'local'
            ])

    def test_run_w_vars(self):
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: echo {{ name }}')

        with self.assertRaises(SystemExit):
            run.main([
                'testwf.jst',
                '--backend',
                'local',
                '--str:name',
                'Philip J. Fry'
            ])

    def test_pipelines(self):
        p = jetstream.Project(new=True)

        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: hostname\n')

        with self.assertRaises(SystemExit):
            pipelines.main([
                'testwf.jst',
                '--backend',
                'local'
            ])

        print(p.workflow().to_yaml(), file=sys.stderr)


    def test_project_tasks(self):
        p = jetstream.Project(new=True)

        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: hostname\n')

        with self.assertRaises(SystemExit):
            pipelines.main([
                'testwf.jst',
                '--backend',
                'local'
            ])

        project.main(['tasks'])
