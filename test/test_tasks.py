from unittest import TestCase
from jetstream.tasks import Task
from jetstream import Workflow


class TaskBasics(TestCase):
    null_task_id_mac = '97d170e1550eee4afc0af065b78cda302a97674c'
    null_task_id = Task().tid

    def test_null_task_id(self):
        """This task will fail if the task identity hashing changes, or
        if there are subtle differences on other platforms. This should not
        affect behavior but it would be nice to understand why it occurs"""
        self.assertEqual(self.null_task_id, self.null_task_id_mac)
    
    def test_task_init(self):
        t = Task(name='taskA')
        self.assertEqual(t.directives()['name'], 'taskA')

    def test_task_init_kw_or_data(self):
        self.assertRaises(
            Exception,
            Task,
            from_data={
            'tid': 'e9e860921f7ca5f33a8f6db90d0cd91924e5a2dc',
            'directives': {'name': 'taskA'}
            },
            name='banana'
        )
    
    def test_task_init2(self):
        """Task rehydration must be from a valid mapping"""
        self.assertRaises(Exception, Task, from_data=42)
        self.assertRaises(Exception, Task, from_data='hello')
        self.assertRaises(Exception, Task, from_data={'no': 'no'})
        self.assertRaises(Exception, Task, from_data={'directives': {'cmd': 'no'}})

    def test_task_init3(self):
        """Attempts to instantiate tasks from "falsey" data should just
        return a null task"""
        self.assertEqual(Task(), Task(from_data={}))
        self.assertEqual(self.null_task_id, Task(from_data={}).tid)

    def test_task_identity(self):
        """Tasks should maintain a single identity across instances"""
        self.assertEqual(Task(), Task())
        self.assertEqual(self.null_task_id, Task().tid)

    def test_null_tasks(self):
        self.assertEqual(Task(), Task(from_data=[]))
        self.assertEqual(self.null_task_id, Task(from_data=[]).tid)
        self.assertEqual(Task(), Task())
        self.assertEqual(self.null_task_id, Task(from_data=0).tid)

    def test_equality(self):
        """Rehydrated tasks should be equal to initted tasks"""
        t1 = Task(from_data={
            'status': 'new',
            'directives': {'name': 'taskA'},
            'state': {}
        })
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
        
    def test_contains(self):
        """Equality handling means we can test for contains"""
        t1 = Task(name='taskA')
        t2 = Task(name='taskB')
        l = [1, 2, 3, t1, 'A', 'B', 'C']
        
        self.assertIn(t1, l)
        self.assertNotIn(t2, l)

    def test_hash_contains(self):
        """Hash method means tasks can be present in containers that require
        hashing. """
        t1 = Task(name='taskA')
        t2 = Task(name='taskB')
        d  = {t1: True, 'something': 42}

        self.assertIn(t1, d)
        self.assertNotIn(t2, d)

    def test_workflow_contains(self):
        t1 = Task(name='taskA')
        t2 = Task(name='taskB')
        wf = Workflow()
        wf.add_task(t1)

        self.assertIn(t1, wf)
        self.assertNotIn(t2, wf)

    def test_is_operator(self):
        t1 = Task(name='taskA')
        t2 = Task(name='taskA')

        self.assertEqual(t1, t2)
        assert t1 is not t2

    def test_task_change_status(self):
        t = Task(name='task')
        self.assertEqual(t.status, 'new')

    def test_task_change_state(self):
        t = Task(name='task')
        t.state['foo'] = 'bar'
        self.assertEqual(t.state['foo'], 'bar')
        t.state.update(foo='baz')
        self.assertEqual(t.state['foo'], 'baz')

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
        t.pending()
        self.assertEqual(t.status, 'pending')
        self.assertTrue(t.is_pending())
        self.assertFalse(t.is_done())

    def test_neg_task_dependency_failure(self):
        t1 = Task(name='task')
        t2 = Task(name='task2', after='task1')
        
        t1.fail()
        self.assertEqual(t2.status, 'new')

    def test_task_serializedeserialze(self):
        t1 = Task(name='task1', cmd='echo hello world', before='foo')
        t2 = Task.deserialize(t1.serialize())

        self.assertEqual(t1, t2)


class TasksInWorkflows(TestCase):
    def test_task_dependency_failure(self):
        """Failing a task should also cause any dependents to fail.
        The task state should indicate which dependency caused the failure"""
        wf = Workflow()
        t1 = wf.new_task(name='task1')
        t2 = wf.new_task(name='task2', after='task1')
        
        t1.fail()

        self.assertEqual(t2.status, 'failed')
        self.assertIn('dependency_failed', t2.state)

    def test_task_dependency_failure_recursive(self):
        """Dependent task failure should propagate to all child tasks"""
        wf = Workflow()

        t1 = wf.new_task(name='task1')
        t2 = wf.new_task(name='task2', after='task1')
        t3 = wf.new_task(name='task3', after='task2')
        t4 = wf.new_task(name='task4', after='task3')
        t1.fail()

        self.assertEqual(t4.status, 'failed')
        self.assertIn('dependency_failed', t4.state)
