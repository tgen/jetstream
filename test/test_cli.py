import os
import sys
import tempfile
import jetstream
from unittest import TestCase
from contextlib import redirect_stdout
from io import StringIO
from jetstream.cli import main as cli_main

jetstream.settings.clear()
jetstream.settings.read(user=False)
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_TEMPLATES = os.path.join(TESTS_DIR, 'templates')
TEST_PIPELINES = os.path.join(TESTS_DIR, 'pipelines')


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
        self.pipeline('foopipe_3')

    def test_foopipe_4(self):
        """pipeline variables are exported and allow files to be accessed"""
        self.pipeline('foopipe_4')


class TestCliModule(TestCase):
    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(TestCliModule, self).setUp()
        self.original_dir = os.getcwd()
        self.temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self.temp_dir.name)

    def tearDown(self):
        os.chdir(self.original_dir)
        self.temp_dir.cleanup()

    def test_init(self):
        """jetstream init should create a project.yaml"""
        args = ['init', '-l', 'silent']
        cli_main(args)
        self.assertTrue(os.path.exists('jetstream/project.yaml'))

    def test_init_w_csv(self):
        """jetstream init -C should auto-detect csv config file"""
        with open('config.csv', 'w') as fp:
            fp.write('foo,bar\nbaz,42\napple,banana')
        args = ['init', '-C', 'config.csv', '--config-file-type', 'csv-nh']
        cli_main(args)
        p = jetstream.Project()
        third_row = p.index['__config_file__'][2]
        self.assertEqual(third_row[0], 'apple')
        self.assertEqual(third_row[1], 'banana')

    def test_init_w_csv_nh(self):
        """jetstream init -C with csv-nh type set"""
        with open('config.csv', 'w') as fp:
            fp.write('foo,bar\nbaz,42\napple,banana')
        args = ['init', '-C', 'config.csv']
        cli_main(args)
        p = jetstream.Project()
        second_row = p.index['__config_file__'][1]
        self.assertEqual(second_row['foo'], 'apple')
        self.assertEqual(second_row['bar'], 'banana')

    def test_init_w_json(self):
        """jetstream init -C with a json config file"""
        with open('config.json', 'w') as fp:
            fp.write('{"foo": "bar", "baz": 42}')
        args = ['init', '-C', 'config.json']
        cli_main(args)
        p = jetstream.Project()
        self.assertEqual(p.index['foo'], 'bar')
        self.assertEqual(p.index['baz'], 42)

    def test_init_w_yaml(self):
        """jetstream init -C with a yaml config file"""
        with open('config.yaml', 'w') as fp:
            fp.write('foo: bar\nbaz: 42')
        args = ['init', '-C', 'config.yaml']
        cli_main(args)
        p = jetstream.Project()
        self.assertEqual(p.index['foo'], 'bar')
        self.assertEqual(p.index['baz'], 42)

    def test_init_w_json_as_yaml(self):
        """specifying the config file type should ignore file extensions"""
        with open('config.yaml', 'w') as fp:
            fp.write('{"foo": "bar", "baz": 42}')
        args = ['init', '-C', 'config.yaml', '--config-file-type', 'json']
        cli_main(args)
        p = jetstream.Project()
        self.assertEqual(p.index['foo'], 'bar')
        self.assertEqual(p.index['baz'], 42)

    def test_reinit(self):
        """running jetstream init again should not affect project id"""
        args = ['init',]
        cli_main(args)
        p = jetstream.Project()

        args = ['init', '-c', 'foo', 'bar']
        cli_main(args)
        p2 = jetstream.Project()

        self.assertEqual(p.info['id'], p2.info['id'])

    def test_run(self):
        """run simple workflow with localbackend"""
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: "true"\n  stdout: /dev/null\n')

        args = ['run', '--backend', 'local', 'testwf.jst',]
        cli_main(args)

    def test_run_w_vars(self):
        """run simple wokrflow with a couple variables passed as args"""
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: echo {{ name }}\n  stdout: /dev/null\n')

        args = [
            'run',
            'testwf.jst',
            '--backend', 'local',
            '-c', 'str:name', 'Philip J. Fry',
            '-c', 'bool:ok', 'true',
            '-c', 'int:number', '42',
            '-c', 'float:number2', '3.14'
        ]

        cli_main(args)

    def test_tasks(self):
        """jetstream tasks cmd should show list of tasks"""
        p = jetstream.init()

        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: "true"\n  stdout: /dev/null\n')

        args = ['run', '--backend', 'local', 'testwf.jst', '--project', p.paths.path]
        cli_main(args)

        args = ['tasks', '--project', p.paths.path]
        with redirect_stdout(StringIO()):
            cli_main(args)

    def test_render(self):
        """jetstream render should just render and print the template"""
        render_test = """{% for i in range(3) %} 
                                                                        
                              ##//#/##/##                            
                            ###//##/###//##                          
                           ###//###//##//###                         
                           ##///##///##///##                         
                           ##///##///###//##                         
                           ##///##///###//##                         
            #/##/#         ##///##///###//##                         
           #/#/##/#        ##///##///###//##                         
           #/#/##/#        ##///##///###//##                         
           #/#/##/#        ##///##///###//##                         
           #/#/##/#        ##///##///###//##                         
           #/#/##/#        ##///##///###//##                         
           #/#/##/#        ##///##///###//##        *#/#/#(/#        
           #/#//##/////////##///##///###//##       /#//#//#/(#       
           #/##////////////##///##///###//##       ##/##//#//#       
           #///##############///##///###//##       ##/##//#//#       
            #################///##///###//##       ##/##//#//#       
                           ##///##///###//##       ##/##//#//#       
                           ##///##///###//##       ##/##//#//#       
                           ##///##///###//##       ##/##//#//#       
                           ##///##///###//##########//##//#//#       
                           ##///##///###//##############//#//#       
                           ##///##///###//##/////////////##//#       
                           ##///##///###//################//##       
                           ##///##///###//##//////////////###      
                           ##///##///###//##                         
                           ##///##///###//##                         
                           ##///##///###//##                         
                           ##///##///###//##                         
                           ##///##///###//##                         
                           ##///##///###//##  
                           ##///##///###//##                         
                           ##///##///###//## 
                           ##///##///###//##  
                           ##///##///###//##                         
                           ##///##///###//##
                            
        {% endfor %}
        """
        with open('testwf.jst', 'w') as fp:
            fp.write(render_test)

        args = ['render', 'testwf.jst']
        cli_main(args)

    def test_build(self):
        """jetstream build should just render and build the template"""
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: "true"\n  stdout: /dev/null\n')

        args = ['build', 'testwf.jst']
        cli_main(args)

    def test_settings(self):
        """jetstream settings should give information about settings"""
        args = ['settings', ]
        cli_main(args)
