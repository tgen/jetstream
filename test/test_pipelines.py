import os
import tempfile
import yaml
import jetstream
from unittest import TestCase


TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PIPELINES = os.path.join(TESTS_DIR, 'pipelines')


class TestCaseWithTempDir(TestCase):
    def setUp(self):
        """ Establishes a new temp dir for the tests to use """
        super(TestCaseWithTempDir, self).setUp()
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()


class PipelineCreation(TestCaseWithTempDir):
    def validate(self, data):
        manifest_path = os.path.join(self.temp_dir.name, 'pipeline.yaml')
        main_path = os.path.join(self.temp_dir.name, 'main.jst')

        with open(manifest_path, 'w') as fp:
            yaml.dump(data, fp)

        with open(main_path, 'w'):
            pass

        jetstream.Pipeline(self.temp_dir.name, validate=True)

    def test_create_pipeline_1(self):
        """Only the name field is required in pipeline info"""
        data = {
            '__pipeline__': {
                'name': 'unittest_pipeline'
            }
        }
        self.validate(data)

    def test_create_pipeline_2(self):
        """Numeric versions should work"""
        data = {
            '__pipeline__': {
                'name': 'unittest_pipeline',
                'version': 0.1
            }
        }
        self.validate(data)

    def test_create_pipeline_3(self):
        """String versions should also work"""
        data = {
            '__pipeline__': {
                'name': 'unittest_pipeline',
                'version': '0.1'
            }
        }
        self.validate(data)

    def test_create_pipeline_4(self):
        """Other data should be fine"""
        data = {
            '__pipeline__': {
                'name': 'unittest_pipeline',
                'version': '0.1'
            },
            'foo': 'bar'
        }
        self.validate(data)


class PipelineLookup(TestCaseWithTempDir):
    def test_find_pipeline_1(self):
        """find_pipelines should find the pipelines"""
        gen = jetstream.pipelines.find_pipelines(TEST_PIPELINES)
        pipelines = list(gen)

        # or
        pipelines = jetstream.pipelines.list_pipelines(TEST_PIPELINES)

    def test_get_pipeline_1(self):
        """Get pipeline should be able to lookup a pipeline by name"""
        jetstream.pipelines.get_pipeline('foo', searchpath=[TEST_PIPELINES,])

    def test_get_pipeline_2(self):
        """Get pipeline shold be able to lookup a pipeline by name and version"""
        jetstream.pipelines.get_pipeline('foo', version=0.1, searchpath=[TEST_PIPELINES,])

    def test_get_pipeline_3(self):
        """Get pipeline should fail to lookup a pipeline we dont have"""
        with self.assertRaises(Exception):
            jetstream.pipelines.get_pipeline('fooz', searchpath=[TEST_PIPELINES,])

    def test_get_pipeline_4(self):
        """Get pipeline should fail to lookup a pipeline version we dont have"""
        with self.assertRaises(Exception):
            jetstream.pipelines.get_pipeline('fooz', version=42, searchpath=[TEST_PIPELINES,])

    def test_get_pipeline_5(self):
        """Get pipeline shold be able to lookup a pipeline by name and version string"""
        jetstream.pipelines.get_pipeline('foo', version='0.1', searchpath=[TEST_PIPELINES, ])