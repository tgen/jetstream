import os
import tempfile
import jetstream
from jetstream.test import TimedTestCase


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
        p = jetstream.project_init()
        self.assertIsInstance(p, jetstream.Project)

    def test_loading_project_data_json(self):
        p = jetstream.project_init()

        test_data = {
            "sampleA": {
                "data": [
                    "path_to_sampleA_data1.txt",
                    "path_to_sampleA_data2.txt",
                    "path_to_sampleA_data3.txt"
                ]
            }
        }

        test_data_path = os.path.join(p.config_path, 'samples.json')

        with open(test_data_path, 'w') as fp:
            jetstream.utils.json_dump(test_data, fp)

        p.reload()
        self.assertEqual(p.config['samples'], test_data)

    def test_project_run(self):
        wf = jetstream.Workflow()
        wf.new_task(name='task', cmd='echo Hello World')

        p = jetstream.project_init()
        rc = p.run(wf)

        self.assertEqual(rc, 0)

