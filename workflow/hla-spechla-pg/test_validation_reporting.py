import csv
import json
from pathlib import Path
import tempfile
import unittest
from score_validation import bootstrap, status, write_tsv


class ValidationReportingTests(unittest.TestCase):
    def test_summary_preserves_columns_when_an_arm_has_no_completed_runs(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'summary.tsv'
            write_tsv(p,[dict(arm='failed',completed=0),dict(arm='valid',completed=1,global_edits=12)])
            with p.open() as stream:
                rows=list(csv.DictReader(stream,delimiter='\t'))
            self.assertEqual(rows[0]['global_edits'],'')
            self.assertEqual(rows[1]['global_edits'],'12')

    def test_unvalidated_marker_is_invalid(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); (p/'COMPLETE').write_text('')
            self.assertEqual(status(p,'donor')[0],'invalid')

    def test_failed_manifest_is_not_missing_or_complete(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); (p/'manifest.json').write_text(json.dumps(dict(status='failed',error='timeout')))
            self.assertEqual(status(p,'donor'),('failed','timeout'))

    def test_bootstrap_donor_difference_and_empty_cohort(self):
        self.assertIsNone(bootstrap([]))
        r=bootstrap([2,2,2])
        self.assertEqual((r['mean'],r['lower_95'],r['upper_95']),(2,2,2))
        self.assertEqual(r['donor_count'],3)
