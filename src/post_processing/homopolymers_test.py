"""
Tests for homopolymers.py.
"""

import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from homopolymers import find_homopolymer_runs

from refactor_config import ORIGINAL_GENOME_RECORD


class TestHomopolymers(unittest.TestCase):

    def test_find_homopolymer_runs(self):
        """Tests finding homopolymer runs.
        """
        SHORT_HOMO_RUN = 'AAAAAAA'
        FILLER_1 = 'TAGCTAGCTA'
        HOMO_RUN_1 = 'CCCCCCCCC'
        FILLER_2 = 'AGTCG'
        HOMO_RUN_2 = 'TTTTTTTTTTTTT'
        HOMO_RUN_3 = 'AAAAAAAAAAA'
        HOMO_RUN_4 = 'GGGGGG'
        FILLER_3 = 'CATCGGCTAGCTAGTAG'
        RAW_SEQ = (SHORT_HOMO_RUN + FILLER_1 + HOMO_RUN_1 + FILLER_2 +
                HOMO_RUN_2 + HOMO_RUN_3 + HOMO_RUN_4 + FILLER_3)
        SEQ = Seq(RAW_SEQ, generic_dna)
        SEQ_RECORD = SeqRecord(SEQ)

        # Test both with and without bounds.
        H_RUNS_NO_BOUNDS = find_homopolymer_runs(SEQ_RECORD)
        H_RUNS_WITH_START_BOUND = find_homopolymer_runs(
                SEQ_RECORD, start_bound=len(SHORT_HOMO_RUN))
        H_RUNS_WITH_BOUNDS = find_homopolymer_runs(
                SEQ_RECORD,
                start_bound=len(SHORT_HOMO_RUN),
                end_bound=len(RAW_SEQ) - 2)

        equivalent_h_run_results = [
                H_RUNS_NO_BOUNDS,
                H_RUNS_WITH_START_BOUND,
                H_RUNS_WITH_BOUNDS]

        for h_runs in equivalent_h_run_results:
            self.assertEqual(4, len(h_runs))

            HOMO_RUN_1_START = len(SHORT_HOMO_RUN) + len(FILLER_1)
            HOMO_RUN_1_END = HOMO_RUN_1_START + len(HOMO_RUN_1)
            self.assertEqual(HOMO_RUN_1, h_runs[0]['group'])
            self.assertEqual((HOMO_RUN_1_START, HOMO_RUN_1_END),
                    h_runs[0]['interval'])

            HOMO_RUN_2_START = HOMO_RUN_1_END + len(FILLER_2)
            HOMO_RUN_2_END = HOMO_RUN_2_START + len(HOMO_RUN_2)
            self.assertEqual(HOMO_RUN_2, h_runs[1]['group'])
            self.assertEqual((HOMO_RUN_2_START, HOMO_RUN_2_END),
                    h_runs[1]['interval'])

            self.assertEqual(HOMO_RUN_3, h_runs[2]['group'])
            self.assertEqual(HOMO_RUN_4, h_runs[3]['group'])

        # Test that bounds actually have effect.
        H_RUNS_WITH_TIGHT_START_BOUND = find_homopolymer_runs(
                SEQ_RECORD,
                start_bound=(len(SHORT_HOMO_RUN) + len(FILLER_1) +
                        len(HOMO_RUN_1)))
        self.assertEqual(3, len(H_RUNS_WITH_TIGHT_START_BOUND))

        H_RUNS_WITH_TIGHT_BOUNDS = find_homopolymer_runs(
                SEQ_RECORD,
                start_bound=(len(SHORT_HOMO_RUN) + len(FILLER_1) +
                        len(HOMO_RUN_1)),
                end_bound=(len(RAW_SEQ) - len(FILLER_3) - len(HOMO_RUN_4))
        )
        self.assertEqual(2, len(H_RUNS_WITH_TIGHT_BOUNDS))


    def test_find_homopolymer_runs_real_genome(self):
        """Test finding homopolymers on the full genome.
        """
        h_run_obj_list = find_homopolymer_runs(ORIGINAL_GENOME_RECORD)

        # Check that the results are ordered by start.
        last = h_run_obj_list[0]
        for h_run_obj in h_run_obj_list[1:]:
            self.assertGreater(h_run_obj['interval'][0], last['interval'][0],
                    "last: %s, current: %s" % (str(last), str(h_run_obj)))
            last = h_run_obj


if __name__ == '__main__':
    unittest.main()
