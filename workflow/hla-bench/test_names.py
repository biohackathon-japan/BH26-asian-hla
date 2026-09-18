import unittest
from names import Nomenclature,truth_slot,pred_slot,compare,split_name
NOM=Nomenclature()
def T(x,gene,level='g_group'):return truth_slot(x,level,NOM,gene)
def P(x,gene,level='g_group'):return pred_slot(x,level,NOM,gene)

class NameTests(unittest.TestCase):
    def test_release_pinned(self):
        self.assertEqual(NOM.release,'IPD-IMGT/HLA 3.65.0')
    def test_split(self):
        self.assertEqual(split_name('HLA-A*02:01:01:02L'),('A',('02','01','01','02')))
        self.assertEqual(split_name('14:01','DRB1'),('DRB1',('14','01')))
        self.assertIsNone(split_name('garbage','A'))
    def test_drb1_1401_1454_same_g_group(self):
        truth=[T('14:01','DRB1'),T('09:01','DRB1')]
        self.assertEqual(compare([P('HLA-DRB1*14:54:01','DRB1'),P('DRB1*09:01:02','DRB1')],truth),(1,1,2))
    def test_different_g_groups_fail(self):
        truth=[T('14:44','DRB1'),T('14:01','DRB1')]
        self.assertEqual(compare([P('14:07','DRB1'),P('14:05','DRB1')],truth)[1],0)
    def test_two_field_level_keeps_1401_1454_distinct(self):
        truth=[T('14:01','DRB1','two_field'),T('09:01','DRB1','two_field')]
        self.assertEqual(compare([P('DRB1*14:54','DRB1','two_field'),P('DRB1*09:01','DRB1','two_field')],truth)[1],0)
    def test_ambiguous_prediction_needs_every_alternative(self):
        truth=[T('02:01','A'),T('11:01','A')]
        self.assertEqual(compare([P('A*02:01,A*24:02','A'),P('A*11:01','A')],truth),(1,0,1))
        self.assertEqual(compare([P('A*02:01:01,A*02:01:02','A'),P('A*11:01','A')],truth),(1,1,2))
    def test_truth_ambiguity_any_member(self):
        truth=[T('32:01:01/24:02','A'),T('11:01','A')]
        self.assertEqual(compare([P('A*24:02','A'),P('A*11:01','A')],truth)[1],1)
    def test_no_call_and_unknown(self):
        truth=[T('02:01','A'),T('11:01','A')]
        self.assertEqual(compare([None,P('A*11:01','A')],truth),(0,0,1))
        self.assertIsNone(P('A*99:99','A'))
        self.assertIsNone(T('99:99','A'))
    def test_renamed_truth_member_translated(self):
        # A*02:01:08 (releases 3.0.0-3.47.0) is now A*02:1040:01, a different two-field type
        self.assertEqual(T('02:01:08','A','two_field'),frozenset(['02:1040']))
        self.assertTrue(T('02:01:08','A'))
    def test_deleted_member_dropped_not_whole_slot(self):
        self.assertEqual(T('47:01:01:01/47:01','B','two_field'),frozenset(['47:01']))
    def test_unordered(self):
        truth=[T('02:01','A'),T('11:01','A')]
        self.assertEqual(compare([P('A*11:01','A'),P('A*02:01','A')],truth),(1,1,2))

if __name__=='__main__':unittest.main()
