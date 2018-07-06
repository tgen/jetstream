import sys
import unittest
import logging
import random
import datetime
import jetstream

log = logging.getLogger(__name__)
RESULTS = {}

class TimedTestCase(unittest.TestCase):
    def setUp(self):
        test_id = self.id()

        if not test_id in RESULTS:
            RESULTS[test_id] = []

        self._trial_start_time = datetime.datetime.now()

    def tearDown(self):
        end_time = datetime.datetime.now()
        elapsed = end_time - self._trial_start_time
        test_id = self.id()

        if self._outcome.success:
            result = 'Pass'
        else:
            result = 'Fail'

        RESULTS[test_id].append({
            'result': result,
            'elapsed': elapsed,
        })


class TaskMethods(TimedTestCase):
    def test_task_init(self):
        jetstream.workflows.Task('task')

    def test_task_change_state(self):
        t = jetstream.workflows.Task('task')

        self.assertEqual(t.state, 0)
        self.assertEqual(t.status, 'new')

    #TODO test task methods


class WorkflowBasics(TimedTestCase):
    def test_workflow_repr(self):
        wf = jetstream.Workflow()

        wf.__repr__()

    def test_workflow_pretty(self):
        wf = jetstream.Workflow()

        wf.add_task('hello_world')

        wf.pretty()

    def test_add_task(self):
        wf = jetstream.Workflow()

        wf.add_task('task')

    def test_add_duplicate_task(self):
        wf = jetstream.Workflow()

        wf.add_task('task')

        self.assertRaises(ValueError, wf.add_task, 'task')

    def test_workflow_len(self):
        """ Test the Workflow.__len__() method"""
        wf = jetstream.Workflow()

        wf.add_task('hello')

        self.assertEqual(len(wf), 1)

    def test_get_task(self):
        wf = jetstream.Workflow()

        wf.add_task('task')
        t = wf.get_task('task')

        self.assertEqual(t.id, 'task')

    def test_add_task_w_cmd(self):
        wf = jetstream.Workflow()
        cmd = 'echo hello world'

        wf.add_task('task', cmd=cmd)
        t = wf.get_task('task')

        self.assertEqual(t.cmd, cmd)

    def test_add_task_w_stdin(self):
        wf = jetstream.Workflow()
        stdin = '#!/bin/bash\necho hello world'

        wf.add_task('task', stdin=stdin)
        t = wf.get_task('task')

        self.assertEqual(t.stdin, stdin.encode('utf8'))

    def test_add_task_w_stdout(self):
        wf = jetstream.Workflow()
        stdout = 'out.txt'

        wf.add_task('task', stdout=stdout)
        t = wf.get_task('task')

        self.assertEqual(t.stdout, stdout)
        self.assertEqual(t.stderr, stdout)

    def test_add_task_w_stderr(self):
        wf = jetstream.Workflow()
        stderr = 'err.txt'

        wf.add_task('task', stdout=stderr)
        t = wf.get_task('task')

        self.assertEqual(t.stderr, stderr)

    def test_add_task_w_stdout_stderr(self):
        wf = jetstream.Workflow()
        stdout = 'out.txt'
        stderr = 'err.txt'

        wf.add_task('task', stdout=stdout, stderr=stderr)
        t = wf.get_task('task')

        self.assertEqual(t.stdout, stdout)
        self.assertEqual(t.stderr, stderr)

    def test_add_task_w_cpus(self):
        wf = jetstream.Workflow()
        cpus = 4

        wf.add_task('task', cpus=cpus)
        t = wf.get_task('task')

        self.assertEqual(t.cpus, cpus)

    def test_add_task_w_mem(self):
        wf = jetstream.Workflow()
        mem = 4

        wf.add_task('task', mem=mem)
        t = wf.get_task('task')

        self.assertEqual(t.mem, mem)

    def test_find_by_id(self):
        wf = jetstream.Workflow()

        wf.add_task('task')

        self.assertEqual(wf.find_by_id('t.*'), {'task',})

    def test_find_by_id_fallback(self):
        wf = jetstream.Workflow()

        self.assertEqual(wf.find_by_id('.*', fallback=None), None)

    def test_neg_find_by_id(self):
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.find_by_id, 't.*')

    def test_find_by_ouput(self):
        wf = jetstream.Workflow()

        wf.add_task('task', output='log.txt')
        matches = wf.find_by_output('log.txt')

        self.assertEqual(matches, {'task',})

    def test_find_by_output_fallback(self):
        wf = jetstream.Workflow()

        self.assertEqual(wf.find_by_output('.*', fallback=None), None)

    def test_neg_find_by_output(self):
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.find_by_output, '.*')

    def test_stress_add_tasks(self):
        wf = jetstream.Workflow()

        for i in range(500):
            wf.add_task(str(i))

        self.assertEqual(len(wf), 500)

    def test_stress_add_tasks_transaction(self):
        wf = jetstream.Workflow()

        with wf:
            for i in range(500):
                wf.add_task(str(i))

        self.assertEqual(len(wf), 500)


class WorkflowDependencies(TimedTestCase):
    def test_neg_add_task_w_input(self):
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.add_task, 'task', input='a')

    def test_add_task_w_output(self):
        wf = jetstream.Workflow()

        wf.add_task('task', output='banana.txt')
        t = wf.get_task('task')

        self.assertEqual(type(t.output), list)

    def test_add_task_w_list_output(self):
        wf = jetstream.Workflow()

        wf.add_task('task', output=['input1.txt', 'input2.txt'])
        t = wf.get_task('task')

        self.assertEqual(type(t.output), list)

    def test_add_task_w_input(self):
        wf = jetstream.Workflow()

        wf.add_task('task1', output='banana.txt')
        wf.add_task('task2', input='banana.txt')
        t = wf.get_task('task2')

        self.assertEqual(type(t.input), list)

    def test_add_task_w_list_input(self):
        wf = jetstream.Workflow()

        wf.add_task('task1', output='log1.txt')
        wf.add_task('task2', output='log2.txt')
        wf.add_task('task3', input=['log1.txt', 'log2.txt'])
        t = wf.get_task('task3')

        self.assertEqual(type(t.input), list)

    def test_neg_add_task_w_after(self):
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.add_task, 'task', after='task1')

    def test_neg_add_task_w_after2(self):
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.add_task, 'task', after='task')

    def test_add_task_w_after(self):
        wf = jetstream.Workflow()

        wf.add_task('task1', )
        t1 = wf.get_task('task1')
        wf.add_task('task2', after='task1')
        deps = list(wf.dependencies('task2'))

        self.assertEqual(deps, [t1,])

    def test_neg_add_task_w_before(self):
        wf = jetstream.Workflow()

        self.assertRaises(ValueError, wf.add_task, 'task', before='task1')

    def test_add_task_w_before(self):
        wf = jetstream.Workflow()

        wf.add_task('task1')
        wf.add_task('task2', before='task1')
        t2 = wf.get_task('task2')
        deps = list(wf.dependencies('task1'))

        self.assertEqual(deps, [t2, ])

    def test_is_ready(self):
        wf = jetstream.Workflow()

        wf.add_task('task1')
        wf.add_task('task2', after='task1')

        self.assertEqual(wf.is_ready('task1'), True)
        self.assertEqual(wf.is_ready('task2'), False)

    def test_neg_composition(self):
        wf = jetstream.Workflow()
        wf.add_task('hello')

        wf2 = jetstream.Workflow()
        wf2.add_task('hello', output='hello')
        wf2.add_task('hello2', input='hello')

        self.assertRaises(ValueError, wf.compose, wf2)

    # TODO Pos composition, compose_all


class WorkflowIteration(TimedTestCase):
    def test_workflow_iter(self):
        wf = jetstream.Workflow()
        status = jetstream.workflows.Task.valid_status

        for i in range(100):
            wf.add_task(str(i), status=random.choice((status)))

        i = iter(wf)

        self.assertEqual(type(i), jetstream.workflows.WorkflowIterator)

    def test_workflow_iter_next1(self):
        wf = jetstream.Workflow()

        wf.add_task('task')
        task = wf.get_task('task')

        i = iter(wf)

        self.assertIs(next(i), task)
        self.assertIs(next(i), None)
        self.assertIs(next(i), None)

        task.complete()

        self.assertRaises(StopIteration, next, i)


def summarize_results():
    for test_id, trials in RESULTS.items():
        print('\nTest:', test_id)

        test_n_trials = len(trials)
        test_trial_times = [t['elapsed'] for t in trials]
        test_total_time = sum(test_trial_times, datetime.timedelta(0, 0))
        test_min_time = min(test_trial_times)
        test_max_time = max(test_trial_times)
        test_mean_time = test_total_time / test_n_trials

        print('Mean:', test_mean_time,
              'Min:', test_min_time,
              'Max:', test_max_time)

        print('-' * 50)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        log.critical('Launching timed test suite')
        n = int(sys.argv[1])

        for i in range(n):
            log.critical('Starting trial {}/{}'.format(i, n))
            unittest.main(verbosity=0, exit=False, argv=sys.argv[1:])

        summarize_results()

    else:
        unittest.main(verbosity=0, exit=True)

else:
    workflow = jetstream.Workflow()
    workflow.add_task('task1', cmd='echo hello world')
    workflow.add_task('task2', cmd='okay', after='task1')
    workflow.add_task('task3', output='task3.results')
    workflow.add_task('task4', output='task4.results')
    workflow.add_task('task5', input='.*results')

