import os
import tempfile
import jetstream
from unittest import TestCase
from contextlib import redirect_stdout
from io import StringIO
from jetstream.cli import main as cli_main

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_TEMPLATES = os.path.join(TESTS_DIR, 'templates')
TEST_VARIABLES = os.path.join(TESTS_DIR, 'templates', 'variables.yaml')


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

    def test_should_pass(self):
        """Valid templates should pass when run"""
        templates_dir = os.path.join(TEST_TEMPLATES, 'should_pass')
        templates = os.listdir(templates_dir)

        for f in templates:
            t = os.path.join(templates_dir, f)

            with self.subTest(msg=t):
                args = ['run', '--backend', 'local', t, '-C', TEST_VARIABLES]
                cli_main(args)


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
        render_test = """{% for i in range(tasks) %} 
                                                                        
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

        args = ['render', 'testwf.jst', '-C', TEST_VARIABLES]
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
