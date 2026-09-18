import importlib.util
import itertools
import json
import re
from pathlib import Path
import tempfile
import unittest

import experiment as E
from phase_linkage import Linkage, read_blast


class TestPhaseRepair(unittest.TestCase):
    def test_complementary_pairs_belong_to_same_phase_edge(self):
        path = Path(__file__).parent / 'review-2026-09-18/installed-map_block2_database.py'
        ns = {'__name__': 'test'}
        exec(E.patch_block_linker(path.read_text(), repair_pairing=True), ns)
        class MockMap:
            def main(self, a, b):
                pair = tuple(int(re.search(r'_hap([12])', p)[1]) for p in (a,b))
                class Result:
                    high_score = {(1,1):14, (2,2):23, (2,1):32, (1,2):41}[pair]
                    support_allele = []
                return Result()
        with tempfile.TemporaryDirectory() as d:
            score = Path(d) / 'scores'
            ns.update(Analyze_map=MockMap, blast_map=lambda x: None, outdir=d, score_file=str(score))
            graph = ns['Construct_Graph'](); graph.fragments=['first','second']; graph.get_edge()
            self.assertEqual(score.read_text().splitlines()[1].split()[2], '37;73')

    def test_maximum_and_ties_independent_of_subject_order(self):
        records = [('low', [90, 100, 90, 100]), ('best1', [100, 100, 100, 100]),
                   ('best2', [100, 100, 100, 100])]
        for order in itertools.permutations(records):
            obj = Linkage(dict(order))
            self.assertEqual(obj.high_score, 20000)
            self.assertEqual([x[0] for x in obj.support_allele], ['best1', 'best2'])
        self.assertEqual(Linkage({}).support_allele, [])

    def test_best_hsp_not_last_hsp(self):
        rows = ['q\tA\t99\t100\n', 'q\tA\t100\t10\n', 'q\tB\t90\t100\n']
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'blast.tsv'
            for order in itertools.permutations(rows):
                p.write_text(''.join(order))
                self.assertEqual(read_blast(None, p), {'A': [99, 100], 'B': [90, 100]})

    def test_patch_installed_source_preserves_interface(self):
        path = Path(__file__).parent / 'review-2026-09-18/installed-map_block2_database.py'
        ns = {'__name__': 'test'}
        exec(E.patch_block_linker(path.read_text()), ns)
        result = ns['Linkage']({'low': [90, 100, 90, 100], 'high': [100, 100, 100, 100]})
        self.assertEqual(result.high_score, 20000)
        self.assertTrue(hasattr(ns['Analyze_map'], 'read_blast'))


class TestRunIntegrity(unittest.TestCase):
    def test_changed_config_or_output_cannot_reuse_completion(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'run'; config = {'arm': 'AE-panel'}
            self.assertTrue(E.start(p, config))
            result = p / 'result'; result.write_text('valid')
            E.finish(p, config, {'result': E.sha(result)})
            self.assertTrue(E.completed(p, config))
            with self.assertRaises(ValueError):
                E.completed(p, {'arm': 'AE-ipd365'})
            result.write_text('corrupt')
            with self.assertRaises(ValueError):
                E.completed(p, config)

    def test_partial_run_never_silently_reused(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'run'
            E.start(p, {'arm': 'A-fixed'})
            with self.assertRaises(FileExistsError):
                E.start(p, {'arm': 'A-fixed'})
            E.failed(p, RuntimeError('nested stage failed'))
            self.assertFalse((p / 'COMPLETE').exists())
            self.assertEqual(json.loads((p / 'manifest.json').read_text())['status'], 'failed')


if __name__ == '__main__':
    unittest.main()
