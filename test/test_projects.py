import os
import tempfile
import jetstream
from unittest import TestCase

jetstream.settings.clear()
jetstream.settings.read(user=False)


class ProjectBasics(TestCase):
    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(ProjectBasics, self).setUp()
        self.original_dir = os.getcwd()
        self.temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.temp_dir.name)

    def tearDown(self):
        os.chdir(self.original_dir)
        self.temp_dir.cleanup()

    def test_project_init(self):
        p = jetstream.init()
        self.assertIsInstance(p, jetstream.Project)

    def test_project_init_other_dir(self):
        p = jetstream.init('banana')
        os.path.exists('banana')
        os.path.exists('banana/jetstream/project.yaml')
        self.assertIsInstance(p, jetstream.Project)

    def test_loading_project_data_json(self):
        test_data = {
            "sampleA": {
                "data": [
                    "path_to_sampleA_data1.txt",
                    "path_to_sampleA_data2.txt",
                    "path_to_sampleA_data3.txt"
                ]
            }
        }
        p = jetstream.init(config=test_data)
        self.assertEqual(p.index, test_data)

    def test_project_run(self):
        wf = jetstream.Workflow()
        wf.new_task(name='task', cmd='echo test_project_run', stdout='/dev/null')
        jetstream.init()
        runner = jetstream.runner.Runner()
        runner.start(wf)


class RunnerBasics(TestCase):
    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(RunnerBasics, self).setUp()
        self._original_dir = os.getcwd()
        self._temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self._temp_dir.name)

    def tearDown(self):
        os.chdir(self._original_dir)
        self._temp_dir.cleanup()

    def test_runner(self):
        runner = jetstream.Runner()
        wf = jetstream.Workflow()
        t = wf.new_task(cmd='hostname', stdout='/dev/null')
        runner.start(workflow=wf)

        self.assertTrue(t.is_complete())
    