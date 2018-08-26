import jetstream
import string
from random import Random
from unittest import TestCase

# Seed the rng because otherwise there's a small chance
# of generating a negative control string that matches
# the id pattern.
random = Random(42)


class TestUlids(TestCase):
    def test_pattern_matches_ids(self):
        """Check that our jetstream run_id pattern actually matches the run
        ids that were generating"""
        for i in range(1000):
            id = jetstream.run_id()
            self.assertTrue(jetstream.run_id_pattern.match(id))

    def test_neg_id_pattern_match(self):
        """Make sure the run id pattern does not match other strings."""
        for i in range(1000):
            n = random.randint(0,100)
            id = ''.join([random.choice(string.printable) for _ in range(n)])
            self.assertFalse(jetstream.run_id_pattern.match(id))

