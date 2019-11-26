import string
import jetstream
import tempfile
from random import Random
from unittest import TestCase

# Seed the rng so these tests are deterministic
random = Random(42)


class UtilsTests(TestCase):
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


class UtilsLoaderParserTests(TestCase):
    def run_parser_tst(self, data, fn, expected):
        res = fn(data)
        self.assertEqual(res, expected)

    def run_loader_tst(self, data, fn, expected):
        """Writes given data to a temp file, loads it with the given function 
        and then compares results to expected"""
        temp = tempfile.NamedTemporaryFile()
        with open(temp.name, 'w') as fp:
            fp.write(data)
        res = fn(temp.name)
        self.assertEqual(res, expected)

    def test_csv(self):
        data = "foo,bar\nbaz,42\napple,banana"
        expected = [
            {'foo': 'baz', 'bar': '42'}, 
            {'foo': 'apple', 'bar': 'banana'}
        ]
        parser = jetstream.utils.parse_csv
        loader = jetstream.utils.load_csv
        self.run_parser_tst(data, parser, expected)
        self.run_loader_tst(data, loader, expected)

    def test_parse_csv_nh(self):
        data = "foo,bar\nbaz,42\napple,banana"
        expected = [
            ['foo', 'bar'],
            ['baz', '42'],
            ['apple', 'banana']
        ]
        parser = jetstream.utils.parse_csv_nh
        loader = jetstream.utils.load_csv_nh
        self.run_parser_tst(data, parser, expected)
        self.run_loader_tst(data, loader, expected)

    def test_parse_tsv(self):
        data = "foo\tbar\nbaz\t42\napple\tbanana"
        expected = [
            {'foo': 'baz', 'bar': '42'}, 
            {'foo': 'apple', 'bar': 'banana'}
        ]
        parser = jetstream.utils.parse_tsv
        loader = jetstream.utils.load_tsv
        self.run_parser_tst(data, parser, expected)
        self.run_loader_tst(data, loader, expected)

    def test_parse_tsv_nh(self):
        data = "foo\tbar\nbaz\t42\napple\tbanana"
        expected = [
            ['foo', 'bar'],
            ['baz', '42'],
            ['apple', 'banana']
        ]
        parser = jetstream.utils.parse_tsv_nh
        loader = jetstream.utils.load_tsv_nh
        self.run_parser_tst(data, parser, expected)
        self.run_loader_tst(data, loader, expected)

    def test_parse_tsv_as_txt(self):
        data = "foo\tbar\nbaz\t42\napple\tbanana"
        expected = [
            'foo\tbar',
            'baz\t42',
            'apple\tbanana'
        ]
        parser = jetstream.utils.parse_txt
        loader = jetstream.utils.load_txt
        self.run_parser_tst(data, parser, expected)
        self.run_loader_tst(data, loader, expected)

    def test_parse_single_column_csv(self):
        data = "foo\nbar\nbaz"
        expected = [{'foo': 'bar'}, {'foo': 'baz'}]
        parser = jetstream.utils.parse_csv
        loader = jetstream.utils.load_csv
        self.run_parser_tst(data, parser, expected)
        self.run_loader_tst(data, loader, expected)

    def test_parse_single_column_csv_nh(self):
        data = "foo\nbar\nbaz"
        expected = [['foo'], ['bar'], ['baz']]
        parser = jetstream.utils.parse_csv_nh
        loader = jetstream.utils.load_csv_nh
        self.run_parser_tst(data, parser, expected)
        self.run_loader_tst(data, loader, expected)

    def test_parse_json(self):
        data = '{"foo": "bar", "baz": 42}'
        expected = {'foo': 'bar', 'baz': 42}
        parser = jetstream.utils.parse_json
        loader = jetstream.utils.load_json
        self.run_parser_tst(data, parser, expected)
        self.run_loader_tst(data, loader, expected)

    def test_parse_yaml(self):
        data = 'foo: bar\nbaz: 42\n'
        expected = {'foo': 'bar', 'baz': 42}
        parser = jetstream.utils.parse_yaml
        loader = jetstream.utils.load_yaml
        self.run_parser_tst(data, parser, expected)
        self.run_loader_tst(data, loader, expected)

