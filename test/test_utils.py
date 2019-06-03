import string
import jetstream
from random import Random
from unittest import TestCase

# Seed the rng so these tests are deterministic
random = Random(42)


class UtilsTest(TestCase):
    def test_dict_lookup_dot_notation(self):
        d = {'status': 'new', 'state': {'slurm': {'something': 42}}}
        v = jetstream.utils.dict_lookup_dot_notation(d, 'status')
        self.assertEqual(v, 'new')

        v = jetstream.utils.dict_lookup_dot_notation(d, 'state')
        self.assertEqual(v, {'slurm': {'something': 42}})

        v = jetstream.utils.dict_lookup_dot_notation(d, 'state.slurm.something')
        self.assertEqual(v, 42)

        self.assertRaises(
            KeyError,
            jetstream.utils.dict_lookup_dot_notation, d, 'foo'
        )

        self.assertRaises(
            KeyError,
            jetstream.utils.dict_lookup_dot_notation, d, 'foo.bar'
        )

        self.assertRaises(
            KeyError,
            jetstream.utils.dict_lookup_dot_notation, d, 'state.foo'
        )



