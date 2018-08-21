import random
import jetstream
from jetstream.tasks import Task
from jetstream.test import TimedTestCase


class WorkflowBasics(TimedTestCase):
    def test_workflow_repr(self):
        wf = jetstream.Workflow()
        wf.__repr__()

    def test_workflow_pretty(self):
        wf = jetstream.Workflow()
        wf.new_task(name='hello_world')
        wf.pretty()

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

        self.assertEqual(t.get('stdin'), stdin)

    def test_add_task_w_stdout(self):
        wf = jetstream.Workflow()
        stdout = 'out.txt'
        t = wf.new_task(name='task', stdout=stdout)

        self.assertEqual(t.get('stdout'), stdout)
        self.assertEqual(t.get('stderr'), None)

    def test_add_task_w_stderr(self):
        wf = jetstream.Workflow()
        stderr = 'err.txt'
        t = wf.new_task(name='task', stderr=stderr)

        self.assertEqual(t.get('stderr'), stderr)

    def test_add_task_w_stdout_stderr(self):
        wf = jetstream.Workflow()
        stdout = 'out.txt'
        stderr = 'err.txt'
        t = wf.new_task(name='task', stdout=stdout, stderr=stderr)

        self.assertEqual(t.get('stdout'), stdout)
        self.assertEqual(t.get('stderr'), stderr)

    def test_add_task_w_cpus(self):
        wf = jetstream.Workflow()
        cpus = 4
        t = wf.new_task(name='task', cpus=cpus)
        
        self.assertEqual(t.get('cpus'), cpus)

    def test_add_task_w_mem(self):
        wf = jetstream.Workflow()
        mem = 4
        t = wf.new_task(name='task', mem=mem)

        self.assertEqual(t.get('mem'), mem)

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


class WorkflowDependencies(TimedTestCase):
    def test_neg_add_task_w_input(self):
        """Task dependencies are linked after every wf.new_task or wf.add_task
        unless context manager is used. Adding a task that requires inputs when
        the no task produces that input as output should raise an error."""
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.new_task, name='task', input='a')

    def test_coerce_sequence_from_str(self):
        wf = jetstream.Workflow()
        t = wf.new_task(name='task', output='banana.txt')
        op = jetstream.utils.coerce_sequence(t.get('output'))      
        
        self.assertIsInstance(op, list)

    def test_coerce_sequence_from_list(self):
        wf = jetstream.Workflow()

        t = wf.new_task(name='task', output=['output1.txt', 'output2.txt'])
        op = jetstream.utils.coerce_sequence(t.get('output'))

        self.assertIsInstance(op, list)

    def test_wf_dependents_1(self):
        wf = jetstream.Workflow()

        t1 = wf.new_task(name='task1', output='banana.txt')
        t2 = wf.new_task(name='task2', input='banana.txt')
        
        deps = set(wf.dependents(t1))
        self.assertEqual(deps, {t2,})

    def test_wf_dependencies_1(self):
        wf = jetstream.Workflow()

        t1 = wf.new_task(name='task1', output='log1.txt')
        t2 = wf.new_task(name='task2', output='log2.txt')
        t3 = wf.new_task(name='task3', input=['log1.txt', 'log2.txt'])
        
        deps = list(wf.dependencies(t3))

        self.assertEqual(deps, [t1, t2])

    def test_neg_add_task_w_after(self):
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.new_task, name='task', after='task1')

    def test_neg_self_dependency(self):
        """Tasks cannot depend on themselves"""
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.new_task, name='task', after='task')

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
    
    # This is no longer a problem because names are not unique
    # def test_neg_composition(self):
    #     wf = jetstream.Workflow()
    #     wf.new_task(name='task1')
    # 
    #     wf2 = jetstream.Workflow()
    #     wf2.new_task(name='task1', output='hello')
    #     wf2.new_task(name='task2', input='hello')
    # 
    #     self.assertRaises(ValueError, wf.compose, wf2)

    # TODO Pos composition, compose_all

    def test_dependent_failure(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='hello')
        t2 = wf.new_task(name='goodbye', after='hello')
        t1.fail()

        self.assertEqual(t2.status, 'failed')


class WorkflowIteration(TimedTestCase):
    def test_workflow_iter(self):
        wf = jetstream.Workflow()
        status = jetstream.workflows.Task.valid_status

        for i in range(100):
            wf.new_task(name=str(i), status=random.choice((status)))

        i = iter(wf)

        self.assertEqual(type(i), jetstream.workflows.WorkflowIterator)

    def test_workflow_iter_next1(self):
        wf = jetstream.Workflow()

        t = wf.new_task(name='task')

        i = iter(wf)

        self.assertIs(next(i), t)
        self.assertIs(next(i), None)
        self.assertIs(next(i), None)

        t.complete()

        self.assertRaises(StopIteration, next, i)
