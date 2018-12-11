import os
import tempfile
import jetstream
from random import Random
from jetstream.tasks import Task
from unittest import TestCase

# Seed the rng because otherwise there's a small chance
# of generating a negative control string that matches
# the id pattern.
random = Random(42)


class WorkflowBasics(TestCase):
    def test_workflow_repr(self):
        wf = jetstream.Workflow()
        wf.__repr__()

    def test_workflow_pretty(self):
        wf = jetstream.Workflow()
        wf.new_task(name='hello_world')
        wf.to_yaml()
        wf.to_json()

    def test_add_task(self):
        wf = jetstream.Workflow()
        wf.new_task(name='task')

    def test_add_duplicate_task(self):
        wf = jetstream.Workflow()
        t1 = Task()
        t2 = Task()    
        wf.add_task(t1)
        
        self.assertRaises(ValueError, wf.add_task, t2)

    def test_workflow_len(self):
        """ Test the Workflow.__len__() method"""
        wf = jetstream.Workflow()
        wf.new_task(name='hello')
        
        self.assertEqual(len(wf), 1)

    def test_get_task(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='task')
        t2 = wf.get_task(t1.tid)

        self.assertIs(t1, t2)

    def test_add_task_w_cmd(self):
        wf = jetstream.Workflow()
        cmd = 'echo hello world'
        t1 = wf.new_task(name='task', cmd=cmd)
        t2 = wf.get_task(t1.tid)
        
        self.assertEqual(t1, t2)
        self.assertIs(t1, t2)

    def test_add_task_w_stdin(self):
        wf = jetstream.Workflow()
        stdin = '#!/bin/bash\necho hello world'
        t = wf.new_task(name='task', stdin=stdin)

        self.assertEqual(t.directives().get('stdin'), stdin)

    def test_add_task_w_stdout(self):
        wf = jetstream.Workflow()
        stdout = 'out.txt'
        t = wf.new_task(name='task', stdout=stdout)

        self.assertEqual(t.directives().get('stdout'), stdout)
        self.assertEqual(t.directives().get('stderr'), None)

    def test_add_task_w_stderr(self):
        wf = jetstream.Workflow()
        stderr = 'err.txt'
        t = wf.new_task(name='task', stderr=stderr)

        self.assertEqual(t.directives().get('stderr'), stderr)

    def test_add_task_w_stdout_stderr(self):
        wf = jetstream.Workflow()
        stdout = 'out.txt'
        stderr = 'err.txt'
        t = wf.new_task(name='task', stdout=stdout, stderr=stderr)

        self.assertEqual(t.directives().get('stdout'), stdout)
        self.assertEqual(t.directives().get('stderr'), stderr)

    def test_add_task_w_cpus(self):
        wf = jetstream.Workflow()
        cpus = 4
        t = wf.new_task(name='task', cpus=cpus)
        
        self.assertEqual(t.directives().get('cpus'), cpus)

    def test_add_task_w_mem(self):
        wf = jetstream.Workflow()
        mem = 4
        t = wf.new_task(name='task', mem=mem)

        self.assertEqual(t.directives().get('mem'), mem)

    def test_find_by_id(self):
        wf = jetstream.Workflow()
        t = wf.new_task(name='task')

        self.assertEqual(wf.find('t.*'), {t.tid,})

    def test_find_by_id_fallback(self):
        wf = jetstream.Workflow()

        self.assertEqual(wf.find('.*', fallback=None), None)

    def test_neg_find_by_id(self):
        wf = jetstream.Workflow()
        self.assertRaises(ValueError, wf.find, 't.*')

    def test_find_by_ouput(self):
        wf = jetstream.Workflow()
        t = wf.new_task(name='task', output='log.txt')
        matches = wf.find_by_output('log.txt')

        self.assertEqual(matches, {t.tid,})

    def test_find_by_output_fallback(self):
        wf = jetstream.Workflow()
        
        self.assertEqual(wf.find_by_output('.*', fallback=None), None)

    def test_neg_find_by_output(self):
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.find_by_output, '.*')

    def test_stress_add_tasks(self):
        wf = jetstream.Workflow()

        for i in range(500):
            wf.new_task(name=str(i))

        self.assertEqual(len(wf), 500)

    def test_stress_add_tasks_transaction(self):
        wf = jetstream.Workflow()

        with wf:
            for i in range(500):
                wf.new_task(name=str(i))

        self.assertEqual(len(wf), 500)

    # TODO Test duplicate task add during context manager


class WorkflowDependencies(TestCase):
    def test_neg_add_task_w_input(self):
        """Task dependencies are linked after every wf.new_task or wf.add_task
        unless context manager is used. Adding a task that requires inputs when
        the no task produces that input as output should raise an error."""
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.new_task, name='task', input='a')

    def test_wf_dependents_1(self):
        wf = jetstream.Workflow()

        t1 = wf.new_task(name='task1', output='banana.txt')
        t2 = wf.new_task(name='task2', input='banana.txt')
        
        deps = set(wf.dependents(t1))
        self.assertEqual(deps, {t2,})

    def test_wf_dependents_2(self):
        wf = jetstream.Workflow()

        t1 = wf.new_task(name='task1', output=['banana.txt', 'banana2.txt'])
        t2 = wf.new_task(name='task2', input=['banana.txt', 'banana2.txt'])

        deps = set(wf.dependents(t1))
        self.assertEqual(deps, {t2,})

    def test_wf_dependencies_1(self):
        wf = jetstream.Workflow()

        t1 = wf.new_task(name='task1', output='log1.txt')
        t2 = wf.new_task(name='task2', output='log2.txt')
        t3 = wf.new_task(name='task3', input=['log1.txt', 'log2.txt'])
        
        deps = set(wf.dependencies(t3))
        self.assertEqual(deps, {t1, t2})

    def test_neg_add_task_w_after(self):
        wf = jetstream.Workflow()
        self.assertRaises(ValueError, wf.new_task, name='task', after='task1')

    def test_neg_self_dependency(self):
        """Tasks cannot depend on themselves. The current task will be removed
        from the match results when searching for dependencies"""
        wf = jetstream.Workflow()
        t = wf.new_task(name='task', after='task')

        self.assertEqual([], list(t.dependencies()))

    def test_add_task_w_after(self):
        wf = jetstream.Workflow()

        t1 = wf.new_task(name='task1', )
        t2 = wf.new_task(name='task2', after='task1')
        
        deps = list(wf.dependencies(t2))

        self.assertEqual(deps, [t1,])

    def test_neg_add_task_w_before(self):
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.new_task, name='task', before='task1')

    def test_add_task_w_before(self):
        wf = jetstream.Workflow()

        t1 = wf.new_task(name='task1')
        t2 = wf.new_task(name='task2', before='task1')

        deps = list(wf.dependencies(t1))

        self.assertEqual(deps, [t2, ])

    def test_is_ready(self):
        wf = jetstream.Workflow()

        t1 = wf.new_task(name='task1')
        t2 = wf.new_task(name='task2', after='task1')
        
        self.assertEqual(wf.is_ready(t1), True)
        self.assertEqual(wf.is_ready(t2), False)

    def test_dependent_failure(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='hello')
        t2 = wf.new_task(name='goodbye', after='hello')
        t1.fail()

        self.assertEqual(t2.status, 'failed')

    def test_workflow_mash(self):
        wf1 = jetstream.random_workflow(n=10)
        wf2 = jetstream.random_workflow(n=20)

        common_task = Task(name='in_common')
        wf1.add_task(common_task)
        wf2.add_task(common_task)

        for t in wf1:
            t.complete()

        wf3 = jetstream.workflows.mash(wf1, wf2)
        self.assertIn(common_task, wf3)

        t = wf3.get_task('in_common')
        self.assertTrue(t.is_complete)



class WorkflowIteration(TestCase):
    def test_random_workflow_n(self):
        wf = jetstream.random_workflow(n=25)

    def test_random_workflow_timout(self):
        wf = jetstream.random_workflow(n=None, timeout=1)

    def test_random_workflow_n_and_timout(self):
        wf = jetstream.random_workflow(250, timeout=1)

    def test_workflow_iter(self):
        wf = jetstream.random_workflow(n=25, timeout=1)
        i = iter(wf)

        self.assertEqual(type(i), jetstream.workflows.Workflow)

    def test_workflow_iter_next1(self):
        wf = jetstream.Workflow()

        t = wf.new_task(name='task')

        i = iter(wf)

        self.assertIs(next(i), t)
        self.assertIs(next(i), None)
        self.assertIs(next(i), None)

        t.complete()

        self.assertRaises(StopIteration, next, i)


class WorkflowSaving(TestCase):
    def setUp(self):
        """ All of these tests take place in the context of a project
        directory. So setUp creates a temp dir and chdir to it. """
        super(WorkflowSaving, self).setUp()
        self._original_dir = os.getcwd()
        self._temp_dir = tempfile.TemporaryDirectory()
        os.chdir(self._temp_dir.name)

    def tearDown(self):
        os.chdir(self._original_dir)
        self._temp_dir.cleanup()

    def test_workflow_save(self):
        wf = jetstream.Workflow()

        with wf:
            for i in range(100):
                wf.new_task(name=str(i), cmd='echo task {}'.format(i))

        jetstream.save_workflow(wf, 'wf.yaml')

    def test_workflow_load(self):
        wf = jetstream.Workflow()
        t1 = jetstream.Task(name='whatever')

        with wf:
            wf.add_task(t1)
            for i in range(100):
                wf.new_task(name=str(i), cmd='echo task {}'.format(i))

        jetstream.save_workflow(wf, 'wf.yaml')

        wf2 = jetstream.load_workflow('wf.yaml')

        self.assertIn(t1, wf2)
