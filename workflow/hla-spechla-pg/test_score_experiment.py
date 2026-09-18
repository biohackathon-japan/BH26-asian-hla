import unittest
from score_experiment import gene_score

class TestGlobalSequenceScore(unittest.TestCase):
    def test_swapped_haplotypes(self):
        self.assertEqual(gene_score(['ACGT','AGGT'],['AGGT','ACGT'])['exact_genes'],2)
    def test_terminal_errors_and_masks_cannot_be_hidden(self):
        r=gene_score(['TTACGT','ANGT'],['ACGT','ACGT'])
        self.assertEqual(r['global_edits'],3)
        self.assertEqual(r['legacy_infix_edits'],1)
        self.assertEqual(r['masked_bases'],1)
        self.assertEqual(r['exact_genes'],0)
