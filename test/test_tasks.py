from unittest import TestCase
from jetstream.tasks import Task
from jetstream import Workflow


class TaskBasics(TestCase):
    null_task_id = 'bf21a9e8fbc5a3846fb05b4fa0859e0917b2202f'
    
    def test_task_init(self):
        t = Task(name='taskA')
        self.assertEqual(t.tid, 'f7d5a5fcfe3b0371be46ca1653978224acd41424')
    
    def test_task_init2(self):
        """Task rehydration must be from a mapping"""
        self.assertRaises(Exception, Task, (42,))
        self.assertRaises(Exception, Task, ('hello',))

    def test_null_tasks_0(self):
        self.assertEqual(Task(), Task())
        self.assertEqual(self.null_task_id, Task().tid)
        
    def test_null_tasks_1(self):
        """"Falsey" data objects should be basically ignored. """
        self.assertEqual(Task(), Task(None))
        self.assertEqual(self.null_task_id, Task(None).tid)
        
    def test_null_tasks_2(self):
        self.assertEqual(Task(), Task({}))
        self.assertEqual(self.null_task_id, Task({}).tid)
        
    def test_null_tasks_3(self):
        self.assertEqual(Task(), Task([]))
        self.assertEqual(self.null_task_id, Task([]).tid)
    
    def test_null_tasks_4(self):
        self.assertEqual(Task(), Task(0))
        self.assertEqual(self.null_task_id, Task(0).tid)
    
    def test_equality(self):
        """Rehydrated tasks should be equal to initted tasks"""
        t1 = Task({'name': 'taskA'})
        t2 = Task(name='taskA')
        self.assertEqual(t1, t2)
        self.assertEqual(t1.tid, t2.tid)
    
    def test_equality2(self):
        """State should not change equality"""
        t1 = Task(name='taskA')
        t2 = Task(name='taskA')
        t1.complete()
        self.assertEqual(t1, t2)
        self.assertEqual(t1.tid, t2.tid)
        
    def test_inclusion(self):
        """Equality handling means we can test for inclusions"""
        t1 = Task(name='taskA')
        t2 = Task(name='taskA')
        l = [t1, 1, 2, 3, 'A', 'B', 'C']
        
        self.assertIn(t2, l)
        
    def test_inclusion_2(self):
        """Inclusion tests can be broken with voodoo, so be careful"""
        class Obj:
            tid = self.null_task_id
            
        x = Obj()
        x.tid = self.null_task_id
        
        y = Task()
        l = [x]
        
        self.assertIn(y, l)

    def test_task_change_state(self):
        t = Task(name='task')
        self.assertEqual(t.status, 'new')

    def test_complete_task(self):
        t = Task(name='task')
        t.complete()
        self.assertEqual(t.status, 'complete')

    def test_fail_task(self):
        t = Task(name='task')
        t.fail()
        self.assertEqual(t.status, 'failed')

    def test_pending_task(self):
        t = Task(name='task')
        t.start()
        self.assertEqual(t.status, 'pending')
        self.assertTrue(t.is_pending())
        self.assertFalse(t.is_done())

    def test_neg_task_dependency_failure(self):
        t1 = Task(name='task')
        t2 = Task(name='task2', after='task1')
        
        t1.fail()
        self.assertEqual(t2.status, 'new')


class TasksInWorkflows(TestCase):
    def test_task_dependency_failure(self):
        wf = Workflow()
        t1 = wf.new_task(name='task1')
        t2 = wf.new_task(name='task2', after='task1')
        
        t1.fail()
        
        self.assertEqual(t2.status, 'failed')

    def test_task_dependency_failure_recursive(self):
        wf = Workflow()

        t1 = wf.new_task(name='task1')
        t2 = wf.new_task(name='task2', after='task1')
        t3 = wf.new_task(name='task3', after='task2')
        t4 = wf.new_task(name='task4', after='task3')
        t1.fail()

        self.assertEqual(t4.status, 'failed')
        self.assertEqual(t4.state['returncode'], 123)
