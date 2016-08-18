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

from vcf_to_genbank import _get_deletion
from vcf_to_genbank import _get_insertion
from vcf_to_genbank import _is_deletion
from vcf_to_genbank import _is_insertion
from vcf_to_genbank import VCFToGenbankMaker


class TestVCFToGenbankMaker(unittest.TestCase):

    def test_insertion__simple(self):
        """Tests insertion.

        Makes sure both the Genbank and the liftover file look correct.
        """
        before_insertion = 'AAAAAAAAAAAAA'
        after_insertion = 'CCCCCCCCCCCCCC'
        raw_seq_str = before_insertion + after_insertion
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        INSERTION_BASE = 'T'

        insertion_data = {
            'position': len(before_insertion),
            'sequence': INSERTION_BASE
        }
        maker.handle_insertion(insertion_data)

        EXPECTED_SEQ = before_insertion + INSERTION_BASE + after_insertion
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))

        # Assert the liftover is correct.
        EXPECTED_LIFTOVER_MAPPING_1 = ((0, len(before_insertion) - 1),
                (0, len(before_insertion) - 1))
        EXPECTED_LIFTOVER_MAPPING_2 = (
                (len(before_insertion), len(raw_seq_str) - 1),
                (len(before_insertion) + 1, len(raw_seq_str)))
        EXPECTED_LIFTOVER = [
                EXPECTED_LIFTOVER_MAPPING_1,
                EXPECTED_LIFTOVER_MAPPING_2]
        self.assertEqual(EXPECTED_LIFTOVER,
                maker.runtime_liftover._interval_mapping)


    def test_insertion__greater_than_one_base(self):
        """Tests insertion with more than one base.
        """
        before_insertion = 'AAAAAAAAAAAAA'
        after_insertion = 'CCCCCCCCCCCCCC'
        raw_seq_str = before_insertion + after_insertion
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        INSERTION_SEQ = 'TT'

        insertion_data = {
            'position': len(before_insertion),
            'sequence': INSERTION_SEQ
        }
        maker.handle_insertion(insertion_data)

        # Assert the resulting sequence is correct.
        EXPECTED_SEQ = before_insertion + INSERTION_SEQ + after_insertion
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))

        # Assert the liftover is correct.
        EXPECTED_LIFTOVER_MAPPING_1 = ((0, len(before_insertion) - 1),
                (0, len(before_insertion) - 1))
        EXPECTED_LIFTOVER_MAPPING_2 = (
                (len(before_insertion), len(raw_seq_str) - 1),
                (len(before_insertion) + 2, len(raw_seq_str) + 1))
        EXPECTED_LIFTOVER = [
                EXPECTED_LIFTOVER_MAPPING_1,
                EXPECTED_LIFTOVER_MAPPING_2]
        self.assertEqual(EXPECTED_LIFTOVER,
                maker.runtime_liftover._interval_mapping)


    def test_insertion__multiple(self):
        """Tests insertion for multiple inserts.

        Makes sure both the Genbank and the liftover file look correct.
        """
        before_insertion_1 = 'AAAAAAAAAAAAA'
        after_insertion_1 = 'CCCCCCCCCCCCCC'
        before_insertion_2 = 'AAAAAAAAAAAAAAAAAA'
        after_insertion_2 = 'CCCCCCCCCCCCCCCC'
        raw_seq_str = (before_insertion_1 + after_insertion_1 +
                before_insertion_2 + after_insertion_2)
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        INSERTION_1_BASE = 'T'
        insertion_1_data = {
            'position': len(before_insertion_1),
            'sequence': INSERTION_1_BASE

        }
        maker.handle_insertion(insertion_1_data)

        INSERTION_2_BASE = 'G'
        insertion_2_data = {
            'position': (len(before_insertion_1) +
                len(after_insertion_1) +
                len(before_insertion_2)),
            'sequence': INSERTION_2_BASE

        }
        maker.handle_insertion(insertion_2_data)

        EXPECTED_SEQ = (before_insertion_1 + INSERTION_1_BASE +
                after_insertion_1 + before_insertion_2 +
                INSERTION_2_BASE + after_insertion_2)
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))


    # TODO: Figure out why this is failing.
    # def test_insertion__overlapping_features(self):
    #     # Example based on intersection of nei and arbB genes in MG1655.
    #     before_overlap = 'GCCCTGGCTGCCAGCA'
    #     after_overlap = 'GCCGACCGCTTCGG'
    #     raw_seq_str = before_overlap + after_overlap
    #     seq = Seq(raw_seq_str, generic_dna)
    #     seq_record = SeqRecord(seq)

    #     feature_1_loc = FeatureLocation(0,
    #             len(before_overlap), strand=1)
    #     feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
    #     seq_record.features.append(feature_1)

    #     feature_2_loc = FeatureLocation(len(before_overlap), len(raw_seq_str),
    #             strand=1)
    #     feature_2 = SeqFeature(feature_2_loc, type='CDS', id='2')
    #     seq_record.features.append(feature_2)

    #     maker = VCFToGenbankMaker(seq_record, None, None)
    #     overlap_replacement = 'TTAA'
    #     maker._update_genome_record_for_variant(
    #             len(before_overlap), '', overlap_replacement)

    #     # Assert the sequence is correct.
    #     EXPECTED_SEQ = before_overlap + overlap_replacement + after_overlap
    #     self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))

    #     # Assert the feature annotations are still correct.
    #     EXPECTED_FEATURE_1_SEQ = before_overlap + overlap_replacement
    #     self.assertEqual(EXPECTED_FEATURE_1_SEQ,
    #             str(feature_1.extract(seq_record.seq)))

    #     EXPECTED_FEATURE_2_SEQ = overlap_replacement + after_overlap
    #     self.assertEqual(EXPECTED_FEATURE_2_SEQ,
    #             str(feature_2.extract(seq_record.seq)))


    def test_deletion__simple(self):
        """Tests a single deletion.
        """
        before_deletion = 'AAAAA'
        deleted_base = 'T'
        after_deletion = 'GGGGGGGGGGG'
        raw_seq_str = before_deletion + deleted_base + after_deletion
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        deletion_data = {
            'interval': (len(before_deletion), len(before_deletion) + 1)
        }
        maker.handle_deletion(deletion_data)

        # Assert the sequence is correct.
        EXPECTED_SEQ = before_deletion + after_deletion
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))

        # Assert the liftover is correct.
        EXPECTED_LIFTOVER = [
                ((0, 4), (0, 4)),
                ((6, 16), (5, 15))
        ]
        self.assertEqual(EXPECTED_LIFTOVER,
                maker.runtime_liftover._interval_mapping)


    def test_deletion__greater_than_one_base__2_bases(self):
        """Tests a deletion of greater than one base.
        """
        before_deletion = 'AAAAA'
        deleted_bases = 'TT'
        after_deletion = 'GGGGGGGGGGG'
        raw_seq_str = before_deletion + deleted_bases + after_deletion
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        deletion_data = {
            'interval': (len(before_deletion),
                    len(before_deletion) + len(deleted_bases))
        }
        maker.handle_deletion(deletion_data)

        # Assert the sequence is correct.
        EXPECTED_SEQ = before_deletion + after_deletion
        self.assertEqual(len(EXPECTED_SEQ), len(seq_record.seq))
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))

        # Assert the liftover is correct.
        EXPECTED_LIFTOVER = [
                ((0, 4), (0, 4)),
                ((7, 17), (5, 15))
        ]
        self.assertEqual(EXPECTED_LIFTOVER,
                maker.runtime_liftover._interval_mapping)


    def test_deletion__greater_than_one_base__3_bases(self):
        """Tests a deletion of 3 bases.
        """
        before_deletion = 'AAAAA'
        deleted_bases = 'TTT'
        after_deletion = 'GGGGGGGGGGG'
        raw_seq_str = before_deletion + deleted_bases + after_deletion
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        deletion_data = {
            'interval': (len(before_deletion),
                    len(before_deletion) + len(deleted_bases))
        }
        maker.handle_deletion(deletion_data)

        # Assert the sequence is correct.
        EXPECTED_SEQ = before_deletion + after_deletion
        self.assertEqual(len(EXPECTED_SEQ), len(seq_record.seq))
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))

        # Assert the liftover is correct.
        EXPECTED_LIFTOVER = [
                ((0, 4), (0, 4)),
                ((8, 18), (5, 15))
        ]
        self.assertEqual(EXPECTED_LIFTOVER,
                maker.runtime_liftover._interval_mapping)


    def test_deletion__multiple_positions(self):
        """Tests multiple deletions.
        """
        before_deletion = 'AAAAA'
        deleted_bases = 'TTT'
        after_deletion = 'GGGGGGGGGGG'
        other_deletion = 'C'
        after_other_deletion = 'TTTTT'
        raw_seq_str = (before_deletion + deleted_bases + after_deletion +
                other_deletion + after_other_deletion)
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        deletion_data = {
            'interval': (len(before_deletion),
                    len(before_deletion) + len(deleted_bases))
        }
        maker.handle_deletion(deletion_data)


        other_deletion_start = (len(before_deletion) + len(deleted_bases) +
                len(after_deletion))
        other_deletion_end = other_deletion_start + len(other_deletion)
        other_deletion_data = {
            'interval': (other_deletion_start, other_deletion_end)
        }
        maker.handle_deletion(other_deletion_data)

        EXPECTED_SEQ = before_deletion + after_deletion + after_other_deletion
        self.assertEqual(len(EXPECTED_SEQ), len(seq_record.seq))
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))


    def test_deletion__overlapping_features(self):
        # Example based on intersection of nei and arbB genes in MG1655.
        before_overlap = 'GCCCTGGCTGCCAGCA'
        overlap = 'CTAG'
        after_overlap = 'GCCGACCGCTTCGG'
        raw_seq_str = before_overlap + overlap + after_overlap
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_loc = FeatureLocation(0,
                len(before_overlap) + len(overlap), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)

        feature_2_loc = FeatureLocation(len(before_overlap), len(raw_seq_str),
                strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id='2')
        seq_record.features.append(feature_2)

        maker = VCFToGenbankMaker(seq_record, None, None)
        maker._update_genome_record_for_variant(
                len(before_overlap), overlap, '')

        # Assert the sequence is correct.
        EXPECTED_SEQ = before_overlap + after_overlap
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))

        # Assert the feature annotations are still correct.
        EXPECTED_FEATURE_1_SEQ = before_overlap
        self.assertEqual(EXPECTED_FEATURE_1_SEQ,
                str(feature_1.extract(seq_record.seq)))

        EXPECTED_FEATURE_2_SEQ = after_overlap
        self.assertEqual(EXPECTED_FEATURE_2_SEQ,
                str(feature_2.extract(seq_record.seq)))


    def test_combo(self):
        """Tests a combination of insertions, deletions, and snps.
        """
        before_deletion = 'AAAAA'
        deleted_bases = 'TTT'
        after_deletion = 'GGGGGGGGGGG'
        after_deletion_with_snp = 'GGGCGGGGGGG'
        insertion = 'AAA'
        after_insertion = 'TTTTT'
        raw_seq_str = (before_deletion + deleted_bases + after_deletion +
                after_insertion)
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        # Make deletion.
        deletion_data = {
            'interval': (len(before_deletion),
                    len(before_deletion) + len(deleted_bases))
        }
        maker.handle_deletion(deletion_data)

        # Make insertion.
        insertion_start = (len(before_deletion) + len(deleted_bases) +
                len(after_deletion))
        insertion_1_data = {
            'position': insertion_start,
            'sequence': insertion

        }
        maker.handle_insertion(insertion_1_data)

        # Make SNP.
        insertion_start = (len(before_deletion) + len(deleted_bases) +
                after_deletion_with_snp.index('C'))
        snp_data = {
            'position': insertion_start,
            'ref': 'G',
            'alt': 'C'
        }
        maker.handle_snp(snp_data)

        EXPECTED_SEQ = (before_deletion + after_deletion_with_snp + insertion +
                after_insertion)
        self.assertEqual(len(EXPECTED_SEQ), len(seq_record.seq))
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))


    def test_combo__reverse_operations(self):
        """Tests a combination of insertions and deletions, similar to above,
        but with operations applied in reverse order.
        """
        before_deletion = 'AAAAA'
        deleted_bases = 'TTT'
        after_deletion = 'GGGGGGGGGGG'
        insertion = 'AAA'
        after_insertion = 'TTTTT'
        raw_seq_str = (before_deletion + deleted_bases + after_deletion +
                after_insertion)
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        insertion_start = (len(before_deletion) + len(deleted_bases) +
                len(after_deletion))
        insertion_1_data = {
            'position': insertion_start,
            'sequence': insertion

        }
        maker.handle_insertion(insertion_1_data)

        deletion_data = {
            'interval': (len(before_deletion),
                    len(before_deletion) + len(deleted_bases))
        }
        maker.handle_deletion(deletion_data)

        EXPECTED_SEQ = (before_deletion + after_deletion + insertion +
                after_insertion)
        self.assertEqual(len(EXPECTED_SEQ), len(seq_record.seq))
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))


    def test_update_genome_record_for_variant(self):
        """Tests handling a vcf record.
        """
        before_variant= 'AAAAA'
        variant = 'TTTT'
        after_variant = 'CCCCCC'
        variant_replacement = 'GGGGGGG'
        raw_seq_str = before_variant + variant + after_variant
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        # Replace the variant
        variant_position = len(before_variant)
        maker._update_genome_record_for_variant(variant_position, variant,
                variant_replacement)

        # Assert the sequence is correct.
        EXPECTED_SEQ = before_variant + variant_replacement + after_variant
        self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))

        # Assert the liftover is correct.
        EXPECTED_LIFTOVER_MAPPING_1 = (
                (0, len(before_variant) - 1),
                (0, len(before_variant) - 1))
        EXPECTED_LIFTOVER_MAPPING_2 = (
                (len(before_variant) + len(variant), len(raw_seq_str) - 1),
                (len(before_variant) + len(variant_replacement),
                        (len(before_variant) + len(variant_replacement) +
                                len(after_variant) - 1)))
        EXPECTED_LIFTOVER = [
                EXPECTED_LIFTOVER_MAPPING_1,
                EXPECTED_LIFTOVER_MAPPING_2]
        self.assertEqual(EXPECTED_LIFTOVER,
                maker.runtime_liftover._interval_mapping)


    # def test_update_genome_record_for_variant__overlapping_features(self):
    #     """Tests handling a record that lands in a region of overlapping
    #     features.
    #     """
    #     # Example based on intersection of nei and arbB genes in MG1655.
    #     before_overlap = 'GCCCTGGCTGCCAGCA'
    #     overlap = 'CTAG'
    #     after_overlap = 'GCCGACCGCTTCGG'
    #     raw_seq_str = before_overlap + overlap + after_overlap
    #     seq = Seq(raw_seq_str, generic_dna)
    #     seq_record = SeqRecord(seq)

    #     feature_1_loc = FeatureLocation(0,
    #             len(before_overlap) + len(overlap), strand=1)
    #     feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
    #     seq_record.features.append(feature_1)

    #     feature_2_loc = FeatureLocation(len(before_overlap), len(raw_seq_str),
    #             strand=1)
    #     feature_2 = SeqFeature(feature_2_loc, type='CDS', id='2')
    #     seq_record.features.append(feature_2)

    #     maker = VCFToGenbankMaker(seq_record, None, None)
    #     overlap_replacement = 'TTAA'
    #     maker._update_genome_record_for_variant(len(before_overlap), overlap,
    #             overlap_replacement)

    #     # Assert the sequence is correct.
    #     EXPECTED_SEQ = before_overlap + overlap_replacement + after_overlap
    #     self.assertEqual(EXPECTED_SEQ, str(seq_record.seq))

    #     # Feature added to represent swap.
    #     # self.assertEqual(3, len(seq_record.features))

    #     # Assert the feature annotations are still correct.
    #     EXPECTED_FEATURE_1_SEQ = before_overlap + overlap_replacement
    #     self.assertEqual(EXPECTED_FEATURE_1_SEQ,
    #             str(feature_1.extract(seq_record.seq)))

    #     EXPECTED_FEATURE_2_SEQ = overlap_replacement + after_overlap
    #     self.assertEqual(EXPECTED_FEATURE_2_SEQ,
    #             str(feature_2.extract(seq_record.seq)))


    def test_convert_target_position_to_source(self):
        """Simple test for converting a position in the target (new) genome
        back to the corresponding position in the original source genome.
        """
        before_insertion = 'AAAAAAAAAAAAA'
        after_insertion = 'CCCCCCCCCCCCCC'
        raw_seq_str = before_insertion + after_insertion
        seq = Seq(raw_seq_str, generic_dna)
        seq_record = SeqRecord(seq)
        maker = VCFToGenbankMaker(seq_record, None, None)

        INSERTION_BASE = 'T'


        INSERTION_POSITION = len(before_insertion)
        insertion_data = {
            'position': INSERTION_POSITION,
            'sequence': INSERTION_BASE
        }
        maker.handle_insertion(insertion_data)

        self.assertEqual(3,
                maker.runtime_liftover.convert_target_position_to_source(3))
        self.assertIsNone(
                maker.runtime_liftover.convert_target_position_to_source(
                        INSERTION_POSITION))
        self.assertEqual(INSERTION_POSITION,
                maker.runtime_liftover.convert_target_position_to_source(
                        INSERTION_POSITION + 1))




class TestHelpers(unittest.TestCase):
    """Tests for various helper methods.
    """

    def test_is_deletion(self):
        """Tests the _is_deletion() method.
        """
        self.assertTrue(_is_deletion('ATTT', 'AT'))
        self.assertTrue(_is_deletion('A', ''))

        self.assertFalse(_is_deletion('', ''))
        self.assertFalse(_is_deletion('', 'T'))
        self.assertFalse(_is_deletion('AGTTT', 'AT'))
        self.assertFalse(_is_deletion('ATG', 'T'))


    def test_get_deletion(self):
        """Tests the _get_deletion() method.
        """
        self.assertEqual('TT', _get_deletion('ATTT', 'AT'))
        self.assertEqual('A', _get_deletion('A', ''))

        self.assertIsNone(_get_deletion('ATTT', 'AGA'))


    def test_is_insertion(self):
        """Tests the _is_deletion() method.
        """
        self.assertTrue(_is_insertion('AT', 'ATTT'))
        self.assertTrue(_is_insertion('', 'A'))
        self.assertTrue(_is_insertion('', 'ATG'))

        self.assertFalse(_is_insertion('', ''))
        self.assertFalse(_is_insertion('T', ''))
        self.assertFalse(_is_insertion('AT', 'AGTTT'))
        self.assertFalse(_is_insertion('T', 'ATG'))


    def test_get_insertion(self):
        """Tests the _is_deletion() method.
        """
        self.assertEqual('TT', _get_insertion('AT', 'ATTT'))
        self.assertEqual('A', _get_insertion('', 'A'))
        self.assertEqual('ATG', _get_insertion('', 'ATG'))

        self.assertIsNone(_get_insertion('T', 'ATG'))


if __name__ == '__main__':
    unittest.main()
