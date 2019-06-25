import jetstream
from unittest import TestCase


class TestKvargs(TestCase):
    def setUp(self):
        self.og = {
            'something': {
                'nested': {
                    'foo': 42,
                    'bar': 24
                }},
            'another': [1,2,3]
        }

    def test_dot_update(self):
        t = self.og.copy()
        t['another'] = 42
        jetstream.utils.dict_update_dot_notation(self.og, 'another', 42)
        self.assertEqual(self.og, t)

    def test_dot_update2(self):
        t = self.og.copy()
        t['something']['nested']['bar'] = 42
        jetstream.utils.dict_update_dot_notation(self.og, 'something.nested.bar', 24)
        self.assertEqual(self.og, t)
