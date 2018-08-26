import os
import tempfile
import logging
import jetstream
from test import TimedTestCase

jetstream.logs.start_logging(level=logging.INFO)


class ProjectBasics(TimedTestCase):
    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(ProjectBasics, self).setUp()
        self._original_dir = os.getcwd()
        self._temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self._temp_dir.name)

    def tearDown(self):
        os.chdir(self._original_dir)
        self._temp_dir.cleanup()

    def test_project_init(self):
        p = jetstream.Project(new=True)
        self.assertIsInstance(p, jetstream.Project)

    def test_loading_project_data_json(self):
        p = jetstream.Project(new=True)

        test_data = {
            "sampleA": {
                "data": [
                    "path_to_sampleA_data1.txt",
                    "path_to_sampleA_data2.txt",
                    "path_to_sampleA_data3.txt"
                ]
            }
        }

        test_data_path = os.path.join(p.config_dir, 'samples.json')

        with open(test_data_path, 'w') as fp:
            jetstream.utils.json_dump(test_data, fp)

        p.reload()
        self.assertEqual(p.config['samples'], test_data)

    def test_project_run(self):
        wf = jetstream.Workflow()
        wf.new_task(name='task', cmd='echo test_project_run ${JETSTREAM_RUN_ID}')
        p = jetstream.Project(new=True)
        runner = jetstream.runner.AsyncRunner()
        rc = runner.start(workflow=wf, project=p)
        self.assertEqual(rc, 0)




class RunnerBasics(TimedTestCase):
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
        runner = jetstream.AsyncRunner()
        wf = jetstream.Workflow()
        wf.new_task(cmd='hostname')

        rc = runner.start(workflow=wf)
        self.assertEqual(rc, 0)
    