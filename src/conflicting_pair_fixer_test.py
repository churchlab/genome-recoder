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
Tests for methods that deal with fixing overlaps.
"""

import copy
import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from biopython_util import get_feature_by_id
from codon_usage_memex import CodonUsageMemex
import conflicting_pair_fixer
from conflicting_pair_fixer import ConflictingPairFixer
from conflicting_pair_fixer import fix_overlap_pair
from conflicting_pair_fixer import _generate_permitted_feature_seq_variants
import conflicting_pair_common
from conflicting_pair_common import does_feature_have_forbidden_codon_in_region
from conflicting_pair_common import identify_conflict_region_of_interest
from refactor_context import RefactorContext


class TestFixOverlaps(unittest.TestCase):

    def _assert_feature_seq(self, feature, seq_record, feature_id_to_seq_map):
        original_seq = feature_id_to_seq_map[feature.id]
        forward_seq = seq_record.seq[
                feature.location.start:feature.location.end]
        self.assertEqual(original_seq, str(forward_seq))

    def test_fix_overlaps_simple(self):
        feature_1_seq = 'ATGTTTGGG'
        feature_2_seq = 'GGGCCCAAAGTA'
        overlap_start_pos = 6
        whole_seq = feature_1_seq[:overlap_start_pos] + feature_2_seq
        overlap_size = len(feature_1_seq) - overlap_start_pos
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_id_to_seq_map = {
            1: feature_1_seq,
            2: feature_2_seq,
        }

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)
        self._assert_feature_seq(feature_1, seq_record, feature_id_to_seq_map)

        feature_2_loc = FeatureLocation(
                overlap_start_pos, overlap_start_pos + len(feature_2_seq),
                strand=-1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)
        self._assert_feature_seq(feature_2, seq_record, feature_id_to_seq_map)

        # Build and use the overlap fixer.
        updated_seq_record = copy.deepcopy(seq_record)
        refactor_context = RefactorContext(updated_seq_record)
        refactor_context.set_forbidden_codon_set(set(['GGG']))
        cpf = ConflictingPairFixer(refactor_context)
        cpf.fix_overlaps()

        EXPECTED_SEQUENCE = feature_1_seq + feature_2_seq

        self.assertEqual(EXPECTED_SEQUENCE, str(updated_seq_record.seq))
        new_feature_1 = get_feature_by_id(updated_seq_record, feature_1.id)
        new_feature_2 = get_feature_by_id(updated_seq_record, feature_2.id)
        self.assertEqual(
                new_feature_1.location.end, new_feature_2.location.start)
        self._assert_feature_seq(
                new_feature_1, updated_seq_record, feature_id_to_seq_map)
        self._assert_feature_seq(
                new_feature_2, updated_seq_record, feature_id_to_seq_map)


    def test_fix_overlaps_third_feature_none_overlap(self):
        feature_1_seq = 'ATGTTTGGG'
        feature_2_seq = 'GGGCCCAAAGTA'
        inter_f2_f3_junk = 'GTAGCTATCTATCTGGTTAAATC'
        feature_3_seq = 'ATGAAACCCTTTGGGTTTCCCAAA'
        overlap_start_pos = 6
        whole_seq = (
                feature_1_seq[:overlap_start_pos] +
                feature_2_seq +
                inter_f2_f3_junk +
                feature_3_seq
        )
        overlap_size = len(feature_1_seq) - overlap_start_pos
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_id_to_seq_map = {
            1: feature_1_seq,
            2: feature_2_seq,
            3: feature_3_seq,
        }

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)
        self._assert_feature_seq(feature_1, seq_record, feature_id_to_seq_map)

        feature_2_loc = FeatureLocation(
                overlap_start_pos, overlap_start_pos + len(feature_2_seq),
                strand=-1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)
        self._assert_feature_seq(feature_2, seq_record, feature_id_to_seq_map)

        feature_3_start = feature_2_loc.end + len(inter_f2_f3_junk)
        feature_3_end = feature_3_start + len(feature_3_seq)
        feature_3_loc = FeatureLocation(
                feature_3_start, feature_3_end, strand=1)
        feature_3 = SeqFeature(feature_3_loc, type='CDS', id=3)
        seq_record.features.append(feature_3)
        self._assert_feature_seq(feature_3, seq_record, feature_id_to_seq_map)

        # Build and use the overlap fixer.
        updated_seq_record = copy.deepcopy(seq_record)
        refactor_context = RefactorContext(updated_seq_record)
        refactor_context.set_forbidden_codon_set(set(['GGG']))
        cpf = ConflictingPairFixer(refactor_context)
        cpf.fix_overlaps()

        EXPECTED_SEQUENCE = (
                feature_1_seq +
                feature_2_seq +
                inter_f2_f3_junk +
                feature_3_seq
        )

        self.assertEqual(EXPECTED_SEQUENCE, str(updated_seq_record.seq))
        for feature_id in feature_id_to_seq_map.keys():
            new_feature = get_feature_by_id(updated_seq_record, feature_id)
            self._assert_feature_seq(
                    new_feature, updated_seq_record, feature_id_to_seq_map)


class TestFixOverlapPair(unittest.TestCase):

    def setUp(self):
        # TODO: Is there a cleaner way to substitute this?
        conflicting_pair_common.RBS_BUFFER_SIZE = 15
        conflicting_pair_fixer.RBS_BUFFER_SIZE = 15
        RBS_BUFFER_SIZE = 15

    def _assert_feature_seq(self, feature, seq_record, feature_id_to_seq_map):
        original_seq = feature_id_to_seq_map[feature.id]
        forward_seq = seq_record.seq[
                feature.location.start:feature.location.end]
        self.assertEqual(original_seq, str(forward_seq))

    def test_fix_overlap_pair_opposing_strands(self):
        """This should simply pull them apart, without adding anything in
        between.
        """
        seq = Seq('ATGTTTGGGCCCAAAGTA', generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_seq = 'ATGTTTGGG'
        feature_2_seq = 'GGGCCCAAAGTA'
        feature_id_to_seq_map = {
            1: feature_1_seq,
            2: feature_2_seq,
        }

        feature_1_loc = FeatureLocation(0, 9, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)
        self._assert_feature_seq(feature_1, seq_record, feature_id_to_seq_map)

        feature_2_loc = FeatureLocation(6, 18, strand=-1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)
        self._assert_feature_seq(feature_2, seq_record, feature_id_to_seq_map)

        updated_seq_record = copy.deepcopy(seq_record)
        is_fix_success = fix_overlap_pair(
                feature_1.id, feature_2.id, updated_seq_record)
        self.assertTrue(is_fix_success)

        self.assertEqual('ATGTTTGGGGGGCCCAAAGTA', str(updated_seq_record.seq))
        new_feature_1 = get_feature_by_id(updated_seq_record, feature_1.id)
        new_feature_2 = get_feature_by_id(updated_seq_record, feature_2.id)
        self.assertEqual(
                new_feature_1.location.end, new_feature_2.location.start)
        self._assert_feature_seq(
                new_feature_1, updated_seq_record, feature_id_to_seq_map)
        self._assert_feature_seq(
                new_feature_2, updated_seq_record, feature_id_to_seq_map)


    def test_fix_overlap_pair_same_direction_forward(self):
        """Account for RBS buffer for right strand.
        """
        overlap_start_pos = 39
        feature_1_seq = 'ATGTTTGGGAAACCCAAACCCGGGTTTAAACCCGGGTTTATGAAAGGG'
        feature_2_seq = feature_1_seq[overlap_start_pos:] + 'CCCAAATTT'
        whole_seq = feature_1_seq[:overlap_start_pos] + feature_2_seq
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_id_to_seq_map = {
            1: feature_1_seq,
            2: feature_2_seq,
        }

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)
        self._assert_feature_seq(feature_1, seq_record, feature_id_to_seq_map)

        feature_2_loc = FeatureLocation(
                overlap_start_pos, overlap_start_pos + len(feature_2_seq),
                strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)
        self._assert_feature_seq(feature_2, seq_record, feature_id_to_seq_map)

        updated_seq_record = copy.deepcopy(seq_record)
        is_fix_success = fix_overlap_pair(
                feature_1.id, feature_2.id, updated_seq_record)
        self.assertTrue(is_fix_success)

        EXPECTED_SEQUENCE = (
                feature_1_seq +
                feature_1_seq[
                        overlap_start_pos - conflicting_pair_common.RBS_BUFFER_SIZE:
                        overlap_start_pos] +
                feature_2_seq
        )
        self.assertEqual(EXPECTED_SEQUENCE, str(updated_seq_record.seq))
        new_feature_1 = get_feature_by_id(updated_seq_record, feature_1.id)
        new_feature_2 = get_feature_by_id(updated_seq_record, feature_2.id)
        self.assertEqual(
                new_feature_1.location.end + conflicting_pair_common.RBS_BUFFER_SIZE,
                new_feature_2.location.start)
        self._assert_feature_seq(
                new_feature_1, updated_seq_record, feature_id_to_seq_map)
        self._assert_feature_seq(
                new_feature_2, updated_seq_record, feature_id_to_seq_map)


    def test_fix_overlap_pair_same_direction_reverse(self):
        """Account for RBS buffer for left strand.
        """
        overlap_start_pos = 12
        feature_1_seq = 'AAACCCGGGTTTCCCAAAGTA'
        feature_2_seq = (feature_1_seq[overlap_start_pos:] +
                'CCCTTTGGGAAACCCAAACCCGGGTTTCCCAAATTTGTA')
        whole_seq = feature_1_seq[:overlap_start_pos] + feature_2_seq
        overlap_seq = feature_1_seq[overlap_start_pos:]
        overlap_size = len(overlap_seq)

        # Sanity check, for visual debug if necessary.
        self.assertEqual(
                'AAACCCGGGTTTCCCAAAGTACCCTTTGGGAAACCCAAACCCGGGTTTCCCAAATTTGTA',
                whole_seq)
        self.assertEqual(
                len(whole_seq),
                len(feature_1_seq) + len(feature_2_seq) - overlap_size)

        # Create the sequence.
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_id_to_seq_map = {
            1: feature_1_seq,
            2: feature_2_seq,
        }

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=-1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)
        self._assert_feature_seq(feature_1, seq_record, feature_id_to_seq_map)

        feature_2_loc = FeatureLocation(
                overlap_start_pos, overlap_start_pos + len(feature_2_seq),
                strand=-1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)
        self._assert_feature_seq(feature_2, seq_record, feature_id_to_seq_map)

        updated_seq_record = copy.deepcopy(seq_record)
        is_fix_success = fix_overlap_pair(
                feature_1.id, feature_2.id, updated_seq_record)
        self.assertTrue(is_fix_success)

        result_seq = str(updated_seq_record.seq)

        EXPECTED_SEQUENCE = (
                feature_1_seq +
                feature_2_seq[overlap_size:overlap_size + conflicting_pair_common.RBS_BUFFER_SIZE] +
                feature_2_seq
        )
        # Sanity check, for visual debug if necessary.
        self.assertEqual(
            'AAACCCGGGTTTCCCAAAGTA'
            'CCCTTTGGGAAACCC'
            'CCCAAAGTACCCTTTGGGAAACCCAAACCCGGGTTTCCCAAATTTGTA',
            EXPECTED_SEQUENCE)

        # Assert correct length. (Simple debug check).
        self.assertEqual(len(EXPECTED_SEQUENCE), len(result_seq))

        # Somewhat redundant piece-wise checks from debugging, but might as well
        # leave them here for the future.

        # Assert the entirety of feature_1 is copied over.
        self.assertEqual(feature_1_seq, result_seq[:len(feature_1_seq)])

        # Assert the RBS extension is properly copied over
        self.assertEqual(
                EXPECTED_SEQUENCE[len(feature_1_seq):
                    len(feature_1_seq) + conflicting_pair_common.RBS_BUFFER_SIZE],
                result_seq[len(feature_1_seq):
                    len(feature_1_seq) + conflicting_pair_common.RBS_BUFFER_SIZE])

        # Assert the entirety of feature_2 is copied over at the end.
        self.assertEqual(
                feature_2_seq,
                result_seq[len(result_seq) - len(feature_2_seq):])

        # Assert the whole thing is correct.
        self.assertEqual(EXPECTED_SEQUENCE, result_seq)

        # Make sure the features are preserved.
        new_feature_1 = get_feature_by_id(updated_seq_record, feature_1.id)
        new_feature_2 = get_feature_by_id(updated_seq_record, feature_2.id)
        self.assertEqual(
                new_feature_1.location.end + conflicting_pair_common.RBS_BUFFER_SIZE,
                new_feature_2.location.start)
        self._assert_feature_seq(
                new_feature_1, updated_seq_record, feature_id_to_seq_map)
        self._assert_feature_seq(
                new_feature_2, updated_seq_record, feature_id_to_seq_map)


    def test_fix_overlap_pair_opposite_directions(self):
        """Account for RBS buffer on both strands.
        """
        overlap_start_pos = 21
        feature_1_seq = 'AAACCCGGGTTTCCCAAACCCATGTTTAAAGGGTTTCCC'
        feature_2_seq = (feature_1_seq[overlap_start_pos:] +
                'CCCTTTGGGAAACCCAAACCCGGGTTTCCCAAATTTAAA')
        whole_seq = feature_1_seq[:overlap_start_pos] + feature_2_seq
        overlap_size = len(feature_1_seq) - overlap_start_pos

        # Sanity check.
        self.assertEqual(
                len(whole_seq),
                len(feature_1_seq) + len(feature_2_seq) - overlap_size)

        # Create the sequence.
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_id_to_seq_map = {
            1: feature_1_seq,
            2: feature_2_seq,
        }

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=-1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)
        self._assert_feature_seq(feature_1, seq_record, feature_id_to_seq_map)

        feature_2_loc = FeatureLocation(
                overlap_start_pos, overlap_start_pos + len(feature_2_seq),
                strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)
        self._assert_feature_seq(feature_2, seq_record, feature_id_to_seq_map)

        updated_seq_record = copy.deepcopy(seq_record)
        is_fix_success = fix_overlap_pair(
                feature_1.id, feature_2.id, updated_seq_record)
        self.assertTrue(is_fix_success)

        result_seq = str(updated_seq_record.seq)

        EXPECTED_SEQUENCE = (
                feature_1_seq +
                feature_2_seq[overlap_size: overlap_size +
                        conflicting_pair_common.RBS_BUFFER_SIZE] +
                feature_1_seq[
                        overlap_start_pos - conflicting_pair_common.RBS_BUFFER_SIZE:
                        overlap_start_pos] +
                feature_2_seq
        )
        self.assertEqual(EXPECTED_SEQUENCE, result_seq)
        new_feature_1 = get_feature_by_id(updated_seq_record, feature_1.id)
        new_feature_2 = get_feature_by_id(updated_seq_record, feature_2.id)
        self.assertEqual(
                new_feature_1.location.end + 2 * conflicting_pair_common.RBS_BUFFER_SIZE,
                new_feature_2.location.start)
        self._assert_feature_seq(
                new_feature_1, updated_seq_record, feature_id_to_seq_map)
        self._assert_feature_seq(
                new_feature_2, updated_seq_record, feature_id_to_seq_map)


class TestIdentifyConflictRegionOfInterest(unittest.TestCase):
    """Tests for method that identifies region of to check for affected
    codons in pairs of features that are overlapping or close enough to affect
    one another.
    """

    def test_identify_conflict_region_of_interest(self):
        feature_1_seq = 'ATGTTTGGG'
        feature_2_seq = 'GGGCCCAAAGTA'
        overlap_start_pos = 6
        whole_seq = feature_1_seq[:overlap_start_pos] + feature_2_seq
        overlap_size = len(feature_1_seq) - overlap_start_pos
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        feature_2_loc = FeatureLocation(
                overlap_start_pos,
                overlap_start_pos + len(feature_2_seq),
                strand=-1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)

        bounds = identify_conflict_region_of_interest(feature_1, feature_2)

        self.assertEqual((6,9), bounds)


    def test_identify_conflict_region_of_interest__none_opposing_no_overlap(self):
        feature_1_seq = 'ATGTTTGGG'
        feature_2_seq = 'GGGCCCAAAGTA'
        whole_seq = feature_1_seq + feature_2_seq
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        feature_2_loc = FeatureLocation(
                len(feature_1_seq),
                len(whole_seq),
                strand=-1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)

        bounds = identify_conflict_region_of_interest(feature_1, feature_2)

        self.assertEqual(None, bounds)


    def test_identify_conflict_region_of_interest__same_dir(self):
        before = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        feature_1_seq = 'ATGTTTGGG'
        gap = 'CCCCCCCCCCC'
        feature_2_seq = 'GGGCCCAAAGTA'
        whole_seq = before + feature_1_seq + gap + feature_2_seq
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_start = len(before)
        feature_1_end = feature_1_start + len(feature_1_seq)
        feature_1_loc = FeatureLocation(
                feature_1_start,
                feature_1_end,
                strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        feature_2_start = feature_1_end + len(gap)
        feature_2_end = feature_2_start + len(feature_2_seq)
        feature_2_loc = FeatureLocation(
                feature_2_start,
                feature_2_end,
                strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)

        bounds = identify_conflict_region_of_interest(feature_1, feature_2)

        EXPECTED_BOUNDS = (feature_2_start - conflicting_pair_common.RBS_BUFFER_SIZE, feature_1_end)
        self.assertEqual(EXPECTED_BOUNDS, bounds)


    def test_identify_conflict_region_of_interest__same_dir_reverse(self):
        before = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        feature_1_seq = 'ATGTTTGGG'
        gap = 'CCCCCCCCCCC'
        feature_2_seq = 'GGGCCCAAAGTA'
        whole_seq = before + feature_1_seq + gap + feature_2_seq
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_start = len(before)
        feature_1_end = feature_1_start + len(feature_1_seq)
        feature_1_loc = FeatureLocation(
                feature_1_start,
                feature_1_end,
                strand=-1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        feature_2_start = feature_1_end + len(gap)
        feature_2_end = feature_2_start + len(feature_2_seq)
        feature_2_loc = FeatureLocation(
                feature_2_start,
                feature_2_end,
                strand=-1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)

        bounds = identify_conflict_region_of_interest(feature_1, feature_2)

        EXPECTED_BOUNDS = (feature_2_start, feature_1_end + conflicting_pair_common.RBS_BUFFER_SIZE)
        self.assertEqual(EXPECTED_BOUNDS, bounds)


    def test_identify_conflict_region_of_interest__head_on(self):
        before = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        feature_1_seq = 'ATGTTTGGG'
        gap = 'CCCCCCCCCCC'
        feature_2_seq = 'GGGCCCAAAGTA'
        whole_seq = before + feature_1_seq + gap + feature_2_seq
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_start = len(before)
        feature_1_end = feature_1_start + len(feature_1_seq)
        feature_1_loc = FeatureLocation(
                feature_1_start,
                feature_1_end,
                strand=-1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        feature_2_start = feature_1_end + len(gap)
        feature_2_end = feature_2_start + len(feature_2_seq)
        feature_2_loc = FeatureLocation(
                feature_2_start,
                feature_2_end,
                strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)

        bounds = identify_conflict_region_of_interest(feature_1, feature_2)

        EXPECTED_BOUNDS = (
                feature_2_start - conflicting_pair_common.RBS_BUFFER_SIZE,
                feature_1_end + conflicting_pair_common.RBS_BUFFER_SIZE
        )
        self.assertEqual(EXPECTED_BOUNDS, bounds)


class TestFixConflictingPairConservative(unittest.TestCase):

    def test_generate_permitted_feature_seq_variants(self):
        CODONS_TO_REMOVE = ['ACC', 'AGG']

        # Simple table for testing.
        AA_TO_CODON_LIST_DICT = {
                'M': {
                        'ATG': {},
                },
                'T': {
                        'ACC': {},
                        'ATT': {},
                        'ACU': {},
                },
                'R': {
                        'AGG': {},
                        'AGA': {},
                        'GGG': {},
                }
        }
        CODON_USAGE_MEMEX = CodonUsageMemex(AA_TO_CODON_LIST_DICT)
        CODON_USAGE_MEMEX.start_codons = ['ATG']

        feature_1_seq = 'ATGACCAGGACU'
        random = 'TTTTCCCTTCGGTT'
        whole_seq = feature_1_seq + random
        seq_record = SeqRecord(whole_seq)

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        # Region covers whole feature.
        region = (0, len(feature_1_seq))
        actual_variants = _generate_permitted_feature_seq_variants(
                feature_1,
                None, # other_feature
                'NA_test', # conflict_type
                seq_record,
                region,
                CODONS_TO_REMOVE,
                CODON_USAGE_MEMEX)
        EXPECTED_VARIANT_SET = set([
                'ATGATTAGAACU', 'ATGATTGGGACU', 'ATGACUAGAACU', 'ATGACUGGGACU',
                'ATGATTAGAATT', 'ATGATTGGGATT', 'ATGACUAGAATT', 'ATGACUGGGATT',
        ])
        self.assertEqual(len(EXPECTED_VARIANT_SET), len(actual_variants))
        self.assertEqual(EXPECTED_VARIANT_SET, set(actual_variants))

        # Region covers partial feature.
        for region_start in [6, 7, 8]:
            region = (region_start, 15)
            actual_variants = _generate_permitted_feature_seq_variants(
                    feature_1,
                    None, # other_feature
                    'NA_test', # conflict_type
                    seq_record,
                    region,
                    CODONS_TO_REMOVE,
                    CODON_USAGE_MEMEX)
            EXPECTED_VARIANT_SET = set([
                'ATGACCAGAACU', 'ATGACCGGGACU', 'ATGACCAGAATT', 'ATGACCGGGATT'
            ])
            self.assertEqual(len(EXPECTED_VARIANT_SET), len(actual_variants))
            self.assertEqual(EXPECTED_VARIANT_SET, set(actual_variants))

        # Region covers inner part of feature.
        for region_start in [6, 7, 8]:
            for region_end in range(region_start + 1, 9):
                region = (region_start, region_end)
                actual_variants = _generate_permitted_feature_seq_variants(
                        feature_1,
                        None, # other_feature
                        'NA_test', # conflict_type
                        seq_record,
                        region,
                        CODONS_TO_REMOVE,
                        CODON_USAGE_MEMEX)
                EXPECTED_VARIANT_SET = set([
                    'ATGACCAGAACU', 'ATGACCGGGACU'
                ])
                self.assertEqual(len(EXPECTED_VARIANT_SET), len(actual_variants))
                self.assertEqual(EXPECTED_VARIANT_SET, set(actual_variants))


class TestFeatureHasForbiddenCodonInRegion(unittest.TestCase):

    def testdoes_feature_have_forbidden_codon_in_region(self):
        FORBIDDEN_CODONS = set(['TAG'])

        before = 'TTTCCCAAATTCC'
        feature_1_seq = 'ATGTAGGGGCCC'
        after = 'TCGATCGATCGT'
        whole_seq = before + feature_1_seq + after
        seq_obj = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq_obj)

        feature_1_loc = FeatureLocation(
                len(before), len(before) + len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        # Make sure the function works for all the possible variations
        # on the region window.
        for region_start in range(len(before), len(before) + 10):
            region = (region_start, region_start + 15)
            if region_start < len(before) + 6:
                assertFn = self.assertTrue
            else:
                assertFn = self.assertFalse
            assertFn(does_feature_have_forbidden_codon_in_region(
                    feature_1, region, seq_record, FORBIDDEN_CODONS))


    def test_reverse_strand(self):
        FORBIDDEN_CODONS = set(['TAG'])

        before = 'TAGTAGTAG'
        feature_1_seq = reverse_complement('ATGTAGGGGCCC')
        after = 'TTTCCCAAATTCC'
        whole_seq = before + feature_1_seq + after
        seq_obj = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq_obj)

        feature_1_loc = FeatureLocation(
                len(before),
                len(before) + len(feature_1_seq),
                strand=-1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)

        # Point-wise checks.
        region = (0, 10)
        self.assertFalse(does_feature_have_forbidden_codon_in_region(
                feature_1, region, seq_record, FORBIDDEN_CODONS))
        region = (3, 9)
        self.assertFalse(does_feature_have_forbidden_codon_in_region(
                feature_1, region, seq_record, FORBIDDEN_CODONS))
        region = (0, 100)
        self.assertTrue(does_feature_have_forbidden_codon_in_region(
                feature_1, region, seq_record, FORBIDDEN_CODONS))
        region = (6, 14)
        self.assertFalse(does_feature_have_forbidden_codon_in_region(
                feature_1, region, seq_record, FORBIDDEN_CODONS))
        region = (6, 15)
        self.assertFalse(does_feature_have_forbidden_codon_in_region(
                feature_1, region, seq_record, FORBIDDEN_CODONS),
                "Relevant part of feature: %s" % (
                    str(seq_record.seq[
                            max(feature_1.location.start, region[0]):
                            min(feature_1.location.end, region[1])])))
        region = (6, 16)
        self.assertTrue(does_feature_have_forbidden_codon_in_region(
                feature_1, region, seq_record, FORBIDDEN_CODONS))

        # Slide an arbitrarily-sized window to check a bunch of cases.
        for region_start in range(0, len(whole_seq)):
            region = (region_start, region_start + 10)
            if region_start >= 6 and region_start <= 17:
                assertFn = self.assertTrue
            else:
                assertFn = self.assertFalse
            assertFn(does_feature_have_forbidden_codon_in_region(
                    feature_1, region, seq_record, FORBIDDEN_CODONS),
                    "Failed for region %s" % str(region))


if __name__ == '__main__':
    unittest.main()
