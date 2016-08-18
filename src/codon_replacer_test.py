# Copyright (C) 2016 The President and Fellows of Harvard College
#
# This file may require additional software or modules to be installed to run
# run properly.
#
# The Genome Recoder software is available under an internal non-commercial
# research and academic use license.  Questions about this software or the
# licensing thereof can be addressed to  Office of Technology Development,
# Harvard University, email: otd@harvard.edu.
#
# @author Gleb Kuznetsov (kuznetsov@g.harvard.edu)

"""
Tests for codon_replacer.py.
"""

import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from codon_replacer import calc_delta_GC_count
from codon_replacer import GraphSearchCodonReplacer
from codon_replacer import SimpleCodonReplacer
from codon_usage_memex import CodonUsageMemex
from feature_profile import GCContentFeatureProfile


class TestSimpleCodonReplacer(unittest.TestCase):
    """Test for the simple codon replacer which just swaps out the first
    synonymous codon available.
    """

    def test_simple_codon_replacer(self):
        CODONS_TO_REMOVE = ['ACC', 'AGG']

        # Simple table for testing.
        AA_TO_CODON_LIST_DICT = {
                'T': {
                        'ACC':{},
                        'ACU': {},
                },
                'R': {
                        'AGG': {},
                        'AGA': {},
                }
        }

        feature_1_seq = 'ATGACCAGG'
        other = 'TTTAAACCCTTT'
        whole_seq = feature_1_seq + other
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        CODON_USAGE_MEMEX = CodonUsageMemex(AA_TO_CODON_LIST_DICT)
        codon_replacer = SimpleCodonReplacer(
                CODONS_TO_REMOVE, CODON_USAGE_MEMEX)

        # Perform replacement.
        replace_result = codon_replacer.replace_codons_in_feature(
                feature_1.id, seq_record)

        EXPECTED_NEW_FEATURE_SEQUENCE = 'ATGACUAGA'
        self.assertEqual(
                EXPECTED_NEW_FEATURE_SEQUENCE,
                str(replace_result['new_feature_seq']))


class TestGraphSearchCodonReplacer(unittest.TestCase):
    """Test for the codon replacer that does a graph search to optimize the
    scores of the profiles.
    """

    def test_graph_search_codon_replacer_gc_content(self):
        CODONS_TO_REMOVE = ['ACC', 'AGG', 'TAG']

        # Simple table for testing.
        AA_TO_CODON_LIST_DICT = {
                'T': {'ACC': {}, 'AAA': {}, 'ACG': {}},
                'R': {'AGG': {}, 'ATT': {}, 'AGC': {}},
                '*': {'TAG': {}, 'TAA': {}, 'AGC': {}},
        }

        prefix = 'CCCAAAGGGCCCAAATTTAAAGGGCCC'
        feature_1_seq = 'ATGACCTAG'
        other = 'TTTAAACCCTTT'
        whole_seq = prefix + feature_1_seq + other
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_start = len(prefix)
        feature_1_end = len(prefix) + len(feature_1_seq)
        feature_1_loc = FeatureLocation(
                feature_1_start, feature_1_end, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        profile_kwargs = {
                'sliding_window_size': 20
        }
        gc_content_profile = GCContentFeatureProfile(
                feature_1.id, seq_record, **profile_kwargs)
        CODON_USAGE_MEMEX = CodonUsageMemex(AA_TO_CODON_LIST_DICT)
        codon_replacer = GraphSearchCodonReplacer(
                CODONS_TO_REMOVE, CODON_USAGE_MEMEX, [gc_content_profile])

        # Perform replacement.
        replace_result = codon_replacer.replace_codons_in_feature(
                feature_1.id, seq_record)

        EXPECTED_NEW_FEATURE_SEQUENCE = 'ATGACGTAA'
        self.assertEqual(
                EXPECTED_NEW_FEATURE_SEQUENCE,
                str(replace_result['new_feature_seq']))


    def test_replacer__start_codons_unchanged(self):
        """Tests that the replacer doesn't swap out forbidden codons in the
        that are serving as a start codon.
        """
        CODONS_TO_REMOVE = ['TTG', 'CTA']

        # Simple table for testing.
        AA_TO_CODON_LIST_DICT = {
                'L': {
                        'TTG': {},
                        'CTT': {},
                        'CTA': {},
                },
                '*': {'TAG': {}, 'TAA': {}, 'AGC': {}},
        }

        CODON_USAGE_MEMEX = CodonUsageMemex(AA_TO_CODON_LIST_DICT)

        feature_1_seq = 'TTG' + 'GCT' + 'AAT' + 'TTG' + 'TAA'
        whole_seq = feature_1_seq
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)

        codon_replacer = GraphSearchCodonReplacer(
                CODONS_TO_REMOVE, CODON_USAGE_MEMEX)

        # Perform replacement.
        replace_result = codon_replacer.replace_codons_in_feature(
                feature_1.id, seq_record)

        # Assert successful fix.
        self.assertTrue(replace_result['is_success'])

        # Assert the new sequence has only the TTG that is not a start codon
        # removed.
        EXPECTED_NEW_FEATURE_SEQUENCE = 'TTG' + 'GCT' + 'AAT' + 'CTT' + 'TAA'
        self.assertEqual(
                EXPECTED_NEW_FEATURE_SEQUENCE,
                str(replace_result['new_feature_seq']))


    def test_replacer__TGA_codon_mid_seq(self):
        """Tests that TGA codon mis-sequence is not replaced,
        as it codes for Selenocysteine.
        """
        CODONS_TO_REMOVE = ['TGA']

        # Simple table for testing.
        AA_TO_CODON_LIST_DICT = {
                '*': {'TGA': {}, 'TAG': {}},
        }

        CODON_USAGE_MEMEX = CodonUsageMemex(AA_TO_CODON_LIST_DICT)

        feature_1_seq = 'TTG' + 'GCT' + 'TGA' + 'TTG' + 'TGA'
        whole_seq = feature_1_seq
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)

        codon_replacer = GraphSearchCodonReplacer(
                CODONS_TO_REMOVE, CODON_USAGE_MEMEX)

        # Perform replacement.
        replace_result = codon_replacer.replace_codons_in_feature(
                feature_1.id, seq_record)

        # Assert successful fix.
        self.assertTrue(replace_result['is_success'])

        # Assert the new sequence has only the TTG that is not a start codon
        # removed.
        EXPECTED_NEW_FEATURE_SEQUENCE = 'TTG' + 'GCT' + 'TGA' + 'TTG' + 'TAG'
        self.assertEqual(
                EXPECTED_NEW_FEATURE_SEQUENCE,
                str(replace_result['new_feature_seq']))


    def test_calc_delta_GC_count(self):
        self.assertEqual(3, calc_delta_GC_count('AAA', 'GGG'))
        self.assertEqual(2, calc_delta_GC_count('AAG', 'GGG'))
        self.assertEqual(1, calc_delta_GC_count('CAC', 'GGG'))
        self.assertEqual(-3, calc_delta_GC_count('GGG', 'AAA'))
        self.assertEqual(-2, calc_delta_GC_count('GCG', 'ACA'))


if __name__ == '__main__':
    unittest.main()
