import jetstream
from unittest import TestCase


class TestKvargs(TestCase):
    def test_kvarg_keys(self):
        x = ['--keyisemptystring:' '--:okay', '--str:better', '--file:csv:best']
        for k in x:
            jetstream.kvargs.parse_key(k)

    def test_nex_kvarg_keys(self):
        x = ['--nope', 'nope', '-nope']
        for k in x:
            self.assertRaises(
                jetstream.kvargs.KvargsError,
                jetstream.kvargs.parse_key,
                k
            )

    def test_kvargs_sequnce(self):
        x = ('--str:hello', 'world', '--int:test', '5')
        d = jetstream.kvargs.parse_kvargs(x)
        self.assertEqual(d, {'hello': 'world', 'test': 5})
