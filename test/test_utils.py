import string
import jetstream
import tempfile
from random import Random
from unittest import TestCase

# Seed the rng so these tests are deterministic
random = Random(42)




test_tsv1 = """foo\tbar
baz\t42
apple\tbanana
"""


class UtilsTests(TestCase):
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
    def parser_test(self, data, fn, expected):
        res = fn(data)
        self.assertEqual(res, expected)

    def loader_test(self, data, fn, expected):
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
        self.parser_test(data, parser, expected)
        self.loader_test(data, loader, expected)

    def test_parse_csv_nh(self):
        data = "foo,bar\nbaz,42\napple,banana"
        expected = [
            ['foo', 'bar'],
            ['baz', '42'],
            ['apple', 'banana']
        ]
        parser = jetstream.utils.parse_csv_nh
        loader = jetstream.utils.load_csv_nh
        self.parser_test(data, parser, expected)
        self.loader_test(data, loader, expected)

    def test_parse_tsv(self):
        data = "foo\tbar\nbaz\t42\napple\tbanana"
        expected = [
            {'foo': 'baz', 'bar': '42'}, 
            {'foo': 'apple', 'bar': 'banana'}
        ]
        parser = jetstream.utils.parse_tsv
        loader = jetstream.utils.load_tsv
        self.parser_test(data, parser, expected)
        self.loader_test(data, loader, expected)

    def test_parse_tsv_nh(self):
        data = "foo\tbar\nbaz\t42\napple\tbanana"
        expected = [
            ['foo', 'bar'],
            ['baz', '42'],
            ['apple', 'banana']
        ]
        parser = jetstream.utils.parse_tsv_nh
        loader = jetstream.utils.load_tsv_nh
        self.parser_test(data, parser, expected)
        self.loader_test(data, loader, expected)

    def test_parse_tsv_as_txt(self):
        data = "foo\tbar\nbaz\t42\napple\tbanana"
        expected = [
            'foo\tbar',
            'baz\t42',
            'apple\tbanana'
        ]
        parser = jetstream.utils.parse_txt
        loader = jetstream.utils.load_txt
        self.parser_test(data, parser, expected)
        self.loader_test(data, loader, expected)

    def test_parse_single_column_csv(self):
        data = "foo\nbar\nbaz"
        expected = [{'foo': 'bar'}, {'foo': 'baz'}]
        parser = jetstream.utils.parse_csv
        loader = jetstream.utils.load_csv
        self.parser_test(data, parser, expected)
        self.loader_test(data, loader, expected)

    def test_parse_single_column_csv_nh(self):
        data = "foo\nbar\nbaz"
        expected = [['foo'], ['bar'], ['baz']]
        parser = jetstream.utils.parse_csv_nh
        loader = jetstream.utils.load_csv_nh
        self.parser_test(data, parser, expected)
        self.loader_test(data, loader, expected)

    def test_parse_json(self):
        data = '{"foo": "bar", "baz": 42}'
        expected = {'foo': 'bar', 'baz': 42}
        parser = jetstream.utils.parse_json
        loader = jetstream.utils.load_json
        self.parser_test(data, parser, expected)
        self.loader_test(data, loader, expected)

    def test_parse_yaml(self):
        data = 'foo: bar\nbaz: 42\n'
        expected = {'foo': 'bar', 'baz': 42}
        parser = jetstream.utils.parse_yaml
        loader = jetstream.utils.load_yaml
        self.parser_test(data, parser, expected)
        self.loader_test(data, loader, expected)

