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
    #longMessage = True

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
        templates_dir = os.path.join(TEST_TEMPLATES, 'should_pass')
        templates = os.listdir(templates_dir)

        for f in templates:
            t = os.path.join(templates_dir, f)

            with self.subTest(msg=t):
                args = ['run', '-l', 'silent', t, '-C', TEST_VARIABLES]
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
        args = ['init', '-l', 'silent']
        cli_main(args)
        self.assertTrue(os.path.exists('jetstream/project.yaml'))

    def test_reinit(self):
        args = ['init', '-l', 'silent']
        cli_main(args)
        p = jetstream.Project()

        args = ['init', '-l', 'silent', '-c', 'foo', 'bar']
        cli_main(args)
        p2 = jetstream.Project()

        self.assertEqual(p.info['id'], p2.info['id'])

    def test_run(self):
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: "true"\n  stdout: /dev/null\n')

        args = ['run', '-l', 'silent', 'testwf.jst',]
        cli_main(args)

    def test_run_w_vars(self):
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: echo {{ name }}\n  stdout: /dev/null\n')

        args = [
            'run',
            'testwf.jst',
            '-l', 'silent',
            '-c', 'str:name', 'Philip J. Fry',
            '-c', 'bool:ok', 'true',
            '-c', 'int:number', '42',
            '-c', 'float:number2', '3.14'
        ]

        cli_main(args)

    def test_tasks(self):
        p = jetstream.init()

        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: "true"\n  stdout: /dev/null\n')

        args = ['run', 'testwf.jst', '-l', 'silent', '--project', p.paths.path]
        cli_main(args)

        args = ['tasks', '--project', p.paths.path]
        with redirect_stdout(StringIO()):
            cli_main(args)

    def test_render(self):
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
        with open('testwf.jst', 'w') as fp:
            fp.write('- cmd: "true"\n  stdout: /dev/null\n')

        args = ['build', 'testwf.jst']
        cli_main(args)

    def test_settings(self):
        args = ['settings', ]
        cli_main(args)
