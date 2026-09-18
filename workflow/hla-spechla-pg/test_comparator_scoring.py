import tempfile
from pathlib import Path
import unittest
from score_comparators import four_field, four_pair, genomic_truth, t1k_calls
from score import compare


class ComparatorScoringTests(unittest.TestCase):
    def test_four_fields_are_not_invented_from_three(self):
        self.assertIsNone(four_field('A*01:01:01'))
        self.assertIsNone(four_pair(['A*01:01:01','A*01:01:01:01'])[0])

    def test_expression_suffix_not_a_numeric_field(self):
        self.assertEqual(four_field('A*01:01:01:02N'),'A*01:01:01:02')

    def test_cds_identity_cannot_supply_genomic_truth(self):
        r=dict(exact_genomic_alleles='',exact_cds_alleles='A*01:01:01:01')
        self.assertIsNone(genomic_truth([r,r]))

    def test_ambiguity_requires_all_predictions_compatible(self):
        truth=(frozenset(['A*01:01:01:01']),frozenset(['A*02:01:01:01']))
        pred=four_pair(['A*01:01:01:01,A*01:01:01:02','A*02:01:01:01'])
        self.assertEqual(compare(pred,truth)[1],0)
        self.assertEqual(compare(four_pair(['A*02:01:01:01','A*01:01:01:01']),truth)[1],1)

    def test_unknown_truth_and_shorter_genomic_label_excluded(self):
        r=dict(exact_genomic_alleles='A*01:01:01:01;A*01:01:01')
        self.assertIsNone(genomic_truth([r,r]))

    def test_t1k_homozygote_and_quality_rule(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'genotype.tsv'
            p.write_text('HLA-A\t1\tHLA-A*01:01:01\t20\t60\t.\t0\t-1\nHLA-B\t2\tHLA-B*07:02:01\t10\t0\tHLA-B*08:01:01\t10\t20\n')
            calls=t1k_calls(p)
            self.assertEqual(calls['A'],['A*01:01:01','A*01:01:01'])
            self.assertEqual(calls['B'],[None,'B*08:01:01'])
