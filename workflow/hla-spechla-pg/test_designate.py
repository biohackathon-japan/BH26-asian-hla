#!/usr/bin/env python3
"""Unit tests for the pangenome designation step (step F)."""
import os,sys,tempfile,unittest
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import designate as D

def levenshtein_hw(q,t):
    """Reference infix edit distance (query q inside target t), for tests without edlib."""
    prev=[0]*(len(t)+1)
    for i,qc in enumerate(q,1):
        cur=[i]+[0]*len(t)
        for j,tc in enumerate(t,1):
            cur[j]=min(prev[j]+1,cur[j-1]+1,prev[j-1]+(qc!=tc))
        prev=cur
    return min(prev)

class TestDesignate(unittest.TestCase):
    def setUp(self):
        self.cat={'S1#1#HLA-A':{'exact_genomic_alleles':'A*11:01:01:01','exact_cds_alleles':'A*11:01:01:01;A*11:01:01:02'},
                  'S2#2#HLA-A':{'exact_genomic_alleles':'','exact_cds_alleles':'A*24:02:01:01;A*24:02:01:02'},
                  'S3#1#HLA-A':{'exact_genomic_alleles':'','exact_cds_alleles':''}}
    def test_gene_body_strips_panel_flanks_only(self):
        seq='F'*D.FLANK+'GENE'+'F'*D.FLANK
        self.assertEqual(D.gene_body('S1#1#HLA-A',seq),'GENE')
        self.assertEqual(D.gene_body('A*01:01:01:01','GENE'),'GENE')
    def test_labels(self):
        self.assertEqual(D.label_for('A*01:01:01:01',self.cat),'A*01:01:01:01')
        self.assertEqual(D.label_for('S1#1#HLA-A',self.cat),'A*11:01:01:01')
        self.assertEqual(D.label_for('S2#2#HLA-A',self.cat),'A*24:02:01:01')  # CDS-level label when no genomic match
        self.assertIsNone(D.label_for('S3#1#HLA-A',self.cat))                 # novel panel sequence
        self.assertEqual(D.two_field('A*24:02:01:01'),'A*24:02');self.assertIsNone(D.two_field(None))
    def test_full_length_filter(self):
        recs=[('A*01:01:01:01','A'*3500),('A*02:01:01:01','A'*3490),('A*03:01:01:01','A'*3505),('A*99:01','A'*1200)]
        self.assertEqual([n for n,_ in D.filter_full_length(recs)],['A*01:01:01:01','A*02:01:01:01','A*03:01:01:01'])
    def test_choose_exact_and_tiebreak(self):
        body='ACGTTGCAAGGCTTACGATCGATCGGATCCTAGCTAGGCTAACGT'*10
        recon='TTTTT'+body+'GGGGG'                                   # flanks around an exact gene body
        alt=body[:200]+('T' if body[200]!='T' else 'A')+body[201:]                                  # one mismatch
        recs=[('A*11:01:01:01',body),('S1#1#HLA-A','N'*D.FLANK+body+'N'*D.FLANK),('A*11:01:99',alt)]
        prior=D.frequency_prior(recs,self.cat)
        win,label,(ee,ge),ties,near=D.choose(recon,recs,self.cat,prior,dist=levenshtein_hw)
        self.assertEqual((ee,ge),(None,0));self.assertEqual(set(ties),{'A*11:01:01:01','S1#1#HLA-A'})
        self.assertEqual(label,'A*11:01:01:01')
        # a novel panel sequence closer than any labelled allele wins and is reported as novel
        novel=body[:100]+('G' if body[100]!='G' else 'C')+body[101:];recs2=recs+[('S3#1#HLA-A','N'*D.FLANK+novel+'N'*D.FLANK)]
        win,label,(ee,ge),ties,near=D.choose('TT'+novel+'TT',recs2,self.cat,prior,dist=levenshtein_hw)
        # the novel record stays the sequence-level winner but the call falls back to the nearest labelled allele
        self.assertEqual((win,ge),('S3#1#HLA-A',0))
        self.assertEqual(label,'A*11:01:01:01');self.assertEqual(near[0],'A*11:01:01:01')
        # with no labelled record in the database at all there is no fallback and no label
        win,label,_,_,near=D.choose('TT'+novel+'TT',[('S3#1#HLA-A','N'*D.FLANK+novel+'N'*D.FLANK)],self.cat,prior,dist=levenshtein_hw)
        self.assertEqual((win,label,near),('S3#1#HLA-A',None,None))
        # tie between a labelled allele and a novel panel record prefers the label
        win,label,_,ties,near=D.choose('TT'+body+'TT',[('A*11:01:01:01',body),('S3#1#HLA-A','N'*D.FLANK+body+'N'*D.FLANK)],self.cat,prior,dist=levenshtein_hw)
        self.assertEqual(win,'A*11:01:01:01')
    def test_exon_stage_precedes_gene_stage(self):
        body='ACGTTGCAAGGCTTACGATCGATCGGATCCTAGCTAGGCTAACGT'*10
        exons={'X*01:01':[(10,60),(200,260)],'X*01:02':[(10,60),(200,260)]}
        # X*01:02 differs from the reconstruction in an intron (3 edits) but not in exons; X*01:01 differs in exon 1 only
        intron_var=body[:120]+'GGG'+body[123:]
        exon_var=body[:30]+('C' if body[30]!='C' else 'G')+body[31:]
        recon='TT'+body+'TT'
        win,label,(ee,ge),ties,near=D.choose(recon,[('X*01:01',exon_var),('X*01:02',intron_var)],{},{},exons=exons,dist=levenshtein_hw)
        self.assertEqual((win,ee),('X*01:02',0));self.assertGreater(ge,0)
        # without exon annotation a record ranks last even when its gene body is closer
        win,label,(ee,ge),ties,near=D.choose(recon,[('X*01:01',exon_var),('X*09:01',body)],{},{},exons={'X*01:01':[(10,60)]},dist=levenshtein_hw)
        self.assertEqual(win,'X*01:01')
    def test_prior_counts_panel_only(self):
        recs=[('S1#1#HLA-A',''),('S2#2#HLA-A',''),('S3#1#HLA-A',''),('A*11:01:01:01','')]
        self.assertEqual(D.frequency_prior(recs,self.cat),{'A*11:01':1,'A*24:02':1})
    def test_edlib_matches_reference_implementation(self):
        try:import edlib
        except ImportError:self.skipTest('edlib not installed locally')
        q='ACGTACGTTGCA';t='GGACGTACCTTGCAGG'
        self.assertEqual(D.edit_distance(q,t),levenshtein_hw(q,t))
    def test_designate_gene_reads_fasta(self):
        body='ACGTTGCAAGGCTTACGATCGATCGGATCCTAGCTAGGCTAACGT'*10
        with tempfile.TemporaryDirectory() as d:
            open(f'{d}/hla.allele.1.HLA_A.fasta','w').write('>HLA_A_0\n'+'CC'+body+'\n')
            recs=[('A*11:01:01:01',body)]
            D.edit_distance=levenshtein_hw
            rows=D.designate_gene(d,'A',recs,self.cat,{})
            self.assertEqual(rows[0]['label'],'A*11:01:01:01');self.assertEqual(rows[0]['gene_edits'],0)
            self.assertEqual(rows[0]['label_source'],'exact')
            self.assertEqual(rows[1]['winner'],'')  # missing hap 2 file -> no call

if __name__=='__main__':unittest.main()
