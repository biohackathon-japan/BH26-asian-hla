import gzip
from pathlib import Path
import tempfile
import unittest
from graph_diagnostic import verify_paths, subset_pairs

class TestGraphGate(unittest.TestCase):
    def test_sequence_mutation_missing_and_leakage_fail(self):
        expected = [('SpecHLA#0#HLA_A','ACGT'), ('donor#1#HLA-A','ATGT')]
        self.assertTrue(verify_paths(expected, [(n+'[0]',s) for n,s in reversed(expected)], set())['passed'])
        self.assertTrue(verify_paths(expected, [(n+'#0',s) for n,s in reversed(expected)], set())['passed'])
        self.assertFalse(verify_paths(expected, expected[:1], set())['passed'])
        self.assertFalse(verify_paths(expected, [expected[0], (expected[1][0],'ACGT')], set())['passed'])
        self.assertFalse(verify_paths(expected, expected, {'donor#1'})['passed'])

    def test_pair_subsetting_rejects_mate_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            for mate in (1,2):
                with gzip.open(d/f'{mate}.gz','wt') as f:
                    f.write(f'@read/{mate}\nAC\n+\nII\n')
            args=[d/'1.gz', d/'2.gz', d/'1.fq', d/'2.fq', 3]
            self.assertEqual(subset_pairs(*args),1)
            with gzip.open(d/'2.gz','wt') as f:
                f.write('@wrong/2\nAC\n+\nII\n')
            with self.assertRaises(ValueError): subset_pairs(*args)
