import random
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from structural_overlay import apply_alleles, project_interval, candidates, long_indels, patch_designator

class TestStructuralOverlay(unittest.TestCase):
    def test_complex_allele_does_not_overwrite_coding_snps(self):
        ref=self.ref[:300]
        altered=('A' if ref[10]!='A' else 'C')
        alt=ref[:10]+altered+ref[11:150]+'N'*60+ref[150:]
        events=long_indels(ref,alt)
        self.assertEqual(len(events),1)
        self.assertEqual(events[0],(149,ref[149],ref[149]+'N'*60))
        result=apply_alleles(ref,ref,events)
        self.assertEqual(result,ref[:150]+'N'*60+ref[150:])
        self.assertEqual(result[10],ref[10])

    def test_full_drb1_window_patch_is_narrow(self):
        self.assertEqual(patch_designator(r'$class\_$j:100-11000 other:100-11000'),r'$class\_$j other:100-11000')
        with self.assertRaises(ValueError): patch_designator('unexpected')

    def test_uncertain_or_heterozygous_graph_alleles_are_not_applied(self):
        class VCF:
            header=SimpleNamespace(samples=['donor'])
            def __enter__(self): return self
            def __exit__(self,*args): pass
            def __iter__(self):
                for gt,gq in [((1,1),0),((0,1),40),((1,1),40)]:
                    yield SimpleNamespace(samples={'donor':dict(GT=gt,GQ=gq,DP=20,AD=(1,19))},
                        alleles=('A','A'+'C'*60),ref='A',alts=('A'+'C'*60,),filter={'PASS'},
                        qual=100,chrom='HLA_DRB1',start=4000,stop=4001,pos=4001)
        with patch.dict('sys.modules',pysam=SimpleNamespace(VariantFile=lambda _: VCF())):
            accepted,audit=candidates('unused','donor','DRB1')
        self.assertEqual(len(accepted),1)
        self.assertEqual([r['decision'] for r in audit],['low_confidence','not_confident_homozygous','accepted'])

    def setUp(self):
        r=random.Random(18)
        self.ref=''.join(r.choice('ACGT') for _ in range(1000))
    def test_insertion_before_replacement_preserved(self):
        seq=self.ref[:100]+'TTTT'+self.ref[100:]
        alt=self.ref[400:430]+'G'*60+self.ref[430:500]
        result=apply_alleles(seq,self.ref,[(400,self.ref[400:500],alt)])
        self.assertEqual(result,seq[:404]+alt+seq[504:])
    def test_deletion_before_replacement_preserved(self):
        seq=self.ref[:100]+self.ref[110:]
        self.assertEqual(project_interval(seq,self.ref,400,500),(390,490))
    def test_ref_mismatch_fails_closed(self):
        with self.assertRaises(ValueError):
            apply_alleles(self.ref,self.ref,[(400,'N'*100,'A'*160)])
    def test_masked_flank_withheld(self):
        seq=self.ref[:370]+'N'*30+self.ref[400:]
        with self.assertRaises(ValueError): project_interval(seq,self.ref,400,500)
    def test_overlap_fails_closed(self):
        with self.assertRaises(ValueError):
            apply_alleles(self.ref,self.ref,[(400,self.ref[400:500],'A'*160),(450,self.ref[450:550],'C'*160)])
