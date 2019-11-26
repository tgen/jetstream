from unittest import TestCase
from jetstream import Workflow, tasks
from jetstream.tasks import Task

jetstream.settings.clear()
jetstream.settings.read(user=False)


class TaskBasics(TestCase):
    def test_init(self):
        Task()

    def test_init2(self):
        Task(name='task')

    def test_init3(self):
        Task(name='task', cmd='echo hello world')

    def test_nonstring_cmd(self):
        self.assertRaises(ValueError, Task, 'taskA', cmd={})

    def test_to_dict(self):
        cmd = 'echo hello world'
        t = Task(name='task', cmd=cmd, foo='bar')
        t.complete()
        d = t.to_dict()
        self.assertEqual(d['name'], 'task')
        self.assertEqual(d['cmd'], cmd)
        self.assertEqual(d['foo'], 'bar')
        self.assertEqual(d['state']['status'], 'complete')

    def test_from_dict(self):
        t1 = tasks.from_dict({
            'name': 'taskA',
            'directives': {'cmd': 'echo hello world'},
            'state': {'status': 'complete'}
        })
        self.assertEqual(t1.name, 'taskA')
        self.assertEqual(t1.status, 'complete')

    def test_neg_deserialize(self):
        """Task rehydration should fail with invalid data"""
        self.assertRaises(Exception, tasks.from_dict, 42)
        self.assertRaises(Exception, tasks.from_dict, None)
        self.assertRaises(Exception, tasks.from_dict, [])
        self.assertRaises(Exception, tasks.from_dict, str)

    def test_identity(self):
        """Equivalent tasks should maintain a single identity across instances
        """
        t1 = Task('foo')
        t2 = Task('foo')

        self.assertEqual(t1, t2)
        self.assertEqual(t1.identity, t2.identity)

        t1 = Task('foo', cmd='echo hello world')
        t2 = Task('foo')

        self.assertEqual(t1, t2)
        self.assertNotEqual(t1.identity, t2.identity)

    def test_equality(self):
        """Rehydrated tasks should be equal to initted tasks"""
        t1 = tasks.from_dict({
            'name': 'taskA',
            'state': {'status': 'new'}
        })
        t2 = Task(name='taskA')
        self.assertEqual(t1, t2)
        self.assertEqual(t1.name, t2.name)
    
    def test_equality2(self):
        """State should not change equality"""
        t1 = Task(name='taskA')
        t2 = Task(name='taskA')
        t1.complete()
        self.assertEqual(t1, t2)
        self.assertEqual(t1.name, t2.name)
        self.assertEqual(t1.identity, t2.identity)
        
    def test_contains(self):
        """Equality handling means we can test for tasks in containers"""
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
        wf.add(t1)
        self.assertIn(t1, wf)
        self.assertNotIn(t2, wf)

    def test_equality_vs_id(self):
        """Tasks can be equal but distinct objects"""
        t1 = Task(name='taskA')
        t2 = Task(name='taskA')
        self.assertEqual(t1, t2)
        self.assertIsNot(t1, t2)
        self.assertEqual(t1.identity, t2.identity)

    def test_change_status(self):
        t = Task(name='task')
        self.assertEqual(t.status, 'new')

    def test_change_state(self):
        t = Task(name='task')
        t.state['foo'] = 'bar'
        self.assertEqual(t.state['foo'], 'bar')
        t.state.update(foo='baz')
        self.assertEqual(t.state['foo'], 'baz')

    def test_complete(self):
        t = Task(name='task')
        t.complete()
        self.assertEqual(t.status, 'complete')
        self.assertTrue(t.is_complete())
        self.assertTrue(t.is_done())
        self.assertFalse(t.is_failed())

    def test_pending(self):
        t = Task()
        t.pending()
        self.assertEqual(t.status, 'pending')
        self.assertTrue(t.is_pending())
        self.assertFalse(t.is_done())

    def test_skipped(self):
        t = Task()
        t.skip(foo='bar')
        self.assertEqual(t.status, 'skipped')
        self.assertTrue(t.is_skipped())
        self.assertTrue(t.is_done())

    def test_fail(self):
        t = Task(name='task')
        t.fail()
        self.assertEqual(t.status, 'failed')

    def test_fail_retry(self):
        t = Task(name='task', retry=3)
        t.fail()
        self.assertEqual(t.status, 'new')
        t.fail()
        self.assertEqual(t.status, 'new')
        t.fail()
        self.assertEqual(t.status, 'new')
        t.fail()
        self.assertEqual(t.status, 'failed')

    def test_pending_task(self):
        t = Task(name='task')
        t.pending()
        self.assertEqual(t.status, 'pending')
        self.assertTrue(t.is_pending())
        self.assertFalse(t.is_done())

    def test_task_serializedeserialze(self):
        t1 = Task(name='task1', cmd='echo hello world', before='foo')
        t2 = tasks.from_dict(t1.to_dict())

        self.assertEqual(t1, t2)
