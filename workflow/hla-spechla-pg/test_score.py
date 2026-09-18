#!/usr/bin/env python3
"""Unit tests for the local scorer."""
import os,sys,tempfile,unittest
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import score as S
from names import Nomenclature,compare

class TestScore(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.nom=Nomenclature()
    def test_read_result_layout(self):
        with tempfile.TemporaryDirectory() as d:
            p=f'{d}/hla.result.txt'
            open(p,'w').write('# version\nSample\tHLA_A_1\tHLA_A_2\tHLA_B_1\tHLA_B_2\nS\tA*01:01:01:01\t-\tB*07:02:01:01\tno_match\n')
            r=S.read_result(p)
            self.assertEqual(r['A'],['A*01:01:01:01',None]);self.assertEqual(r['B'],['B*07:02:01:01',None]);self.assertEqual(r['C'],[None,None])
    def test_three_field(self):
        self.assertEqual(S.three_field('A*01:01:01:01'),'A*01:01:01');self.assertEqual(S.three_field('DRB1*15:02:02'),'DRB1*15:02:02')
        self.assertIsNone(S.three_field('A*01:01'))
    def test_truth_slots_and_compare(self):
        recs=[{'exact_cds_alleles':'A*01:01:01:01;A*01:01:01:02N'},{'exact_cds_alleles':'A*02:01:01:01;A*02:01:01:02'}]
        ts=S.truth_slots(recs,'two_field',self.nom,'A')
        self.assertEqual(ts,(frozenset({'01:01'}),frozenset({'02:01'})))
        ts3=S.truth_slots(recs,'three_field',self.nom,'A')
        self.assertEqual(ts3,(frozenset({'A*01:01:01'}),frozenset({'A*02:01:01'})))
        # unordered pair, no-call is an error
        self.assertEqual(compare(S.pred_pair(['A*02:01:01:01','A*01:01:01:01'],'two_field',self.nom,'A'),ts)[1],1)
        self.assertEqual(compare(S.pred_pair(['A*02:01:01:01',None],'two_field',self.nom,'A'),ts)[1],0)
        self.assertEqual(compare(S.pred_pair(['A*02:01:01:01','A*01:01:02'],'three_field',self.nom,'A'),ts3)[1],0)
        self.assertEqual(compare(S.pred_pair(['A*02:01:01:05','A*01:01:01:99'],'three_field',self.nom,'A'),ts3)[1],1)
        self.assertIsNone(S.truth_slots([{'exact_cds_alleles':''},recs[1]],'two_field',self.nom,'A'))
    def test_g_group_level_accepts_G_names(self):
        recs=[{'exact_cds_alleles':'A*01:01:01:01'},{'exact_cds_alleles':'A*01:01:01:01'}]
        ts=S.truth_slots(recs,'g_group',self.nom,'A')
        self.assertEqual(compare(S.pred_pair(['A*01:01:01G','A*01:01:01G'],'g_group',self.nom,'A'),ts)[1],1)
    def test_sequence_score(self):
        body='ACGTTGCAAGGCTTACGATCGATCGGATCCTAGCTAGGCTAACGT'*20
        alt=body[:300]+('A' if body[300]!='A' else 'C')+body[301:]
        with tempfile.TemporaryDirectory() as d:
            open(f'{d}/hla.allele.1.HLA_A.fasta','w').write('>HLA_A_0\nTTTT'+alt+'GGGG\n')
            open(f'{d}/hla.allele.2.HLA_A.fasta','w').write('>HLA_A_1\nTTTT'+body+'GGGG\n')
            recs=[{'name':'X#1#HLA-A'},{'name':'X#2#HLA-A'}]
            panel={'X#1#HLA-A':'N'*S.FLANK+body+'N'*S.FLANK,'X#2#HLA-A':'N'*S.FLANK+alt+'N'*S.FLANK}
            from pathlib import Path
            e1,e2,exact,nmask=S.sequence_score(Path(d),'A',recs,panel)
            self.assertEqual((e1,e2,exact,nmask),(0,0,2,0))
            os.remove(f'{d}/hla.allele.2.HLA_A.fasta')
            e1,e2,exact,nmask=S.sequence_score(Path(d),'A',recs,panel)
            self.assertEqual(exact,1);self.assertEqual(sorted([e1,e2])[0],0);self.assertEqual(max(e1,e2),len(body))

if __name__=='__main__':unittest.main()
