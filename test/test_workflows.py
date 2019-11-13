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
    def test_init(self):
        jetstream.Workflow()

    def test_len(self):
        wf = jetstream.Workflow()
        self.assertEqual(len(wf), 0)

    def test_add_task(self):
        t = jetstream.Task()
        wf = jetstream.Workflow()
        wf.add(t)

    def test_pop_task(self):
        t = jetstream.Task()
        wf = jetstream.Workflow(tasks=[t,])
        t2 = wf.pop(t.name)
        self.assertIs(t, t2)
        self.assertEqual(len(wf), 0)

    def test_new_task(self):
        wf = jetstream.Workflow()
        wf.new_task(cmd='hello world')
        self.assertEqual(len(wf), 1)

    def test_add_duplicate_task(self):
        wf = jetstream.Workflow()
        t1 = Task()
        t2 = Task()    
        wf.add(t1)
        self.assertRaises(ValueError, wf.add, t2)

    def test_indexing(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task()
        t2 = wf[t1.name]
        t3 = wf[t1.name]
        self.assertIs(t1, t2, t3)

    def test_add_task_w_cmd(self):
        wf = jetstream.Workflow()
        cmd = 'echo hello world'
        t1 = wf.new_task(name='task', cmd=cmd)
        t2 = wf['task']
        self.assertEqual(t1, t2)
        self.assertIs(t1, t2)

    def test_add_task_w_stdin(self):
        wf = jetstream.Workflow()
        stdin = '/path/to/stdin.txt'
        task = wf.new_task(name='task', stdin=stdin)
        self.assertEqual(task.directives.get('stdin'), stdin)

    def test_add_task_w_stdout(self):
        wf = jetstream.Workflow()
        stdout = 'out.txt'
        task = wf.new_task(name='task', stdout=stdout)
        self.assertEqual(task.directives.get('stdout'), stdout)

    def test_add_task_w_stderr(self):
        wf = jetstream.Workflow()
        stderr = 'err.txt'
        task = wf.new_task(name='task', stderr=stderr)
        self.assertEqual(task.directives.get('stderr'), stderr)

    def test_add_task_w_stdout_stderr(self):
        wf = jetstream.Workflow()
        stdout = 'out.txt'
        stderr = 'err.txt'
        task = wf.new_task(name='task', stdout=stdout, stderr=stderr)
        self.assertEqual(task.directives.get('stdout'), stdout)
        self.assertEqual(task.directives.get('stderr'), stderr)

    def test_add_task_w_cpus(self):
        wf = jetstream.Workflow()
        cpus = 4
        task = wf.new_task(name='task', cpus=cpus)
        self.assertEqual(task.directives.get('cpus'), cpus)

    def test_add_task_w_mem(self):
        wf = jetstream.Workflow()
        mem = 4
        task = wf.new_task(name='task', mem=mem)
        self.assertEqual(task.directives.get('mem'), mem)

    def test_find_by_id(self):
        wf = jetstream.Workflow()
        t = wf.new_task(name='task')
        results = wf.find('t.*')
        self.assertIn(t, results)

    def test_find_by_id_fallback(self):
        wf = jetstream.Workflow()
        self.assertEqual(wf.find('.*', fallback=None), None)

    def test_neg_find_by_id(self):
        wf = jetstream.Workflow()
        self.assertRaises(ValueError, wf.find, 't.*')

    def test_stress_add_tasks(self):
        wf = jetstream.Workflow()

        for i in range(500):
            wf.new_task(name=f'task_{i}')

        self.assertEqual(len(wf), 500)

    def test_stress_add_tasks_transaction(self):
        wf = jetstream.Workflow()

        for i in range(500):
            wf.new_task(name=f'task{i}')

        self.assertEqual(len(wf), 500)


class WorkflowDependencies(TestCase):
    def test_successors_1(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='task1', output='banana.txt')
        t2 = wf.new_task(name='task2', input='banana.txt')
        deps = set(wf.graph.successors(t1))
        self.assertEqual(deps, {t2,})

    def test_successors_2(self):
        """Test adding input/output directives as lists"""
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='task1', output=['banana.txt', 'banana2.txt'])
        t2 = wf.new_task(name='task2', input=['banana.txt', 'banana2.txt'])
        deps = set(wf.graph.successors(t1))
        self.assertEqual(deps, {t2,})

    def test_predecessors_1(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='task1', output='log1.txt')
        t2 = wf.new_task(name='task2', output='log2.txt')
        t3 = wf.new_task(name='task3', input=['log1.txt', 'log2.txt'])
        deps = set(wf.graph.predecessors(t3))
        self.assertEqual(deps, {t1, t2})

    def test_neg_self_dependency(self):
        """Tasks cannot depend on themselves. Adding edges where from==to
        does nothing. """
        wf = jetstream.Workflow()
        t = wf.new_task(name='task', after='task')
        deps = set(wf.graph.predecessors(t))
        self.assertEqual(deps, set())

    def test_add_task_w_after(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='task1', )
        t2 = wf.new_task(name='task2', after='task1')
        deps = set(wf.graph.predecessors(t2))
        self.assertEqual(deps, {t1,})

    def test_add_task_w_before(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='task1')
        t2 = wf.new_task(name='task2', before='task1')
        deps = set(wf.graph.predecessors(t1))
        self.assertEqual(deps, {t2,})

    def test_is_ready(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='task1')
        t2 = wf.new_task(name='task2', after='task1')
        self.assertEqual(wf.graph.is_ready(t1), True)
        self.assertEqual(wf.graph.is_ready(t2), False)

    def test_dependent_failure(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='hello')
        t2 = wf.new_task(name='goodbye', after='hello')
        wf.graph.skip_descendants(t1)
        self.assertTrue(t2.is_done())
        self.assertTrue(t2.is_failed())
        self.assertTrue(t2.is_skipped())
        self.assertFalse(t2.is_complete())

    def test_reset_parents(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='taskA')
        t2 = wf.new_task(name='taskB', after='taskA', reset='parents')
        for t in wf:
            t.complete()
        self.assertEqual(t1.status, 'complete')
        wf.reset_task(t2)
        self.assertEqual(t1.status, 'new')

    def test_reset_task_parents_recursive(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='taskA')
        t2 = wf.new_task(name='taskB', after='taskA', reset='parents')
        t3 = wf.new_task(name='taskC', after='taskB', reset='parents')
        for t in wf:
            t.complete()
        self.assertEqual(t1.status, 'complete')
        wf.reset_task(t3)
        self.assertEqual(t1.status, 'new')
    

    def test_reset_task_name(self):
        wf = jetstream.Workflow()
        t1 = wf.new_task(name='taskA')
        t2 = wf.new_task(name='taskB', after='taskA')
        t3 = wf.new_task(name='taskC', after='taskB', reset='taskA')
        for t in wf:
            t.complete()
        self.assertEqual(t1.status, 'complete')
        wf.reset_task(t3)
        self.assertEqual(t1.status, 'new')
    

    def test_workflow_mash(self):
        wf1 = jetstream.random_workflow(n=10)
        wf2 = jetstream.random_workflow(n=20)

        common_task = Task(name='in_common')
        wf1.add(common_task)
        wf2.add(common_task)

        for t in wf1:
            t.complete()

        wf3 = jetstream.workflows.mash(wf1, wf2)
        self.assertIn(common_task, wf3)

        t = wf3['in_common']
        self.assertTrue(t.is_complete())



class WorkflowIteration(TestCase):
    def test_random_workflow_n(self):
        jetstream.random_workflow(n=25)

    def test_random_workflow_timout(self):
        jetstream.random_workflow(n=None, timeout=1)

    def test_random_workflow_n_and_timout(self):
        jetstream.random_workflow(250, timeout=1)

    def test_iter(self):
        wf = jetstream.random_workflow(n=25, timeout=1)
        iter(wf)

    def test_graph(self):
        wf = jetstream.random_workflow(n=25, timeout=1)
        wf.reload_graph()
        wf.graph

    def test_graph_iter(self):
        """While tasks are pending but workflow is not complete, the
        WorkflowGraphIterator will return None. When all tasks are complete
        it should raise StopIteration"""
        wf = jetstream.Workflow()

        t = wf.new_task(name='task')
        i = iter(wf.graph)
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

    def test_save(self):
        wf = jetstream.Workflow()
        for i in range(100):
            wf.new_task(name=f'task{i}', cmd='echo task {}'.format(i))
        jetstream.save_workflow(wf, 'wf.pickle')
        self.assertTrue(os.path.exists('wf.pickle'))

    def test_load(self):
        wf = jetstream.Workflow()

        # Add some tasks
        for i in range(100):
            wf.new_task(name=f'task{i}', cmd='echo task {}'.format(i))

        # Modify some state
        wf['task50'].complete()
        wf['task42'].skip(foo='bar')

        jetstream.save_workflow(wf, 'wf.pickle')
        wf2 = jetstream.load_workflow('wf.pickle')

        self.assertEqual(len(wf), 100)
        self.assertIn('task1', wf2)
        self.assertTrue(wf2['task50'].is_complete())
        self.assertEqual(wf2['task42'].state.get('foo'), 'bar')

