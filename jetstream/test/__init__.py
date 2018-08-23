import sys
import unittest
import logging
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
        log.info('Launching timed test suite')
        n = int(sys.argv[1])

        for i in range(n):
            log.info('Starting trial {}/{}'.format(i, n))
            unittest.main(verbosity=0, exit=False, argv=sys.argv[1:])

        summarize_results()

    else:
        log.info('Running unittests')
        unittest.main(verbosity=0, exit=True)

else:
    workflow = jetstream.Workflow()
    workflow.new_task(name='task1', cmd='echo hello world')
    workflow.new_task(name='task2', cmd='okay', after='task1')
    workflow.new_task(name='task3', output='task3.results')
    workflow.new_task(name='task4', output='task4.results')
    workflow.new_task(name='task5', input='.*results')

