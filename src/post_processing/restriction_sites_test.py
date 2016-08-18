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
Test for the restriction_sites module.
"""

import unittest

from Bio import Restriction
from Bio.Alphabet import generic_dna
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from biopython_util import translate_custom
from feature_profile import CodonRarityFeatureProfile
from feature_profile import GCContentFeatureProfile
from feature_profile import SecondaryStructureFeatureProfile
from refactor_context import RefactorContext
from restriction_sites import _remove_site_in_coding_feature
from restriction_sites import _remove_site_not_in_coding_feature
from restriction_sites import _remove_site_with_partial_feature_overlap
from restriction_sites import find_restriction_site_occurrences
from restriction_sites import remove_enzyme_site


class TestRemoveRestrictionSites(unittest.TestCase):
    """Tests related to removing restriction sites.
    """

    def test_remove_site__coding(self):
        """Test removing a single site that falls in a coding feature.
        """
        pass


    # def test_remove_enzyme_site__non_coding(self):
    #     """Test removing all occurrences of a restriction enzyme site in a
    #     genome.
    #     """
    #     RESTRICTION_ENZYME = Restriction.BsmBI
    #     SITE_SEQ = RESTRICTION_ENZYME.site
    #     BEFORE = 'CCCCCCCCCCCCCCCCCCCCCCCC'
    #     MID = 'TTTTTTTTTTTT'
    #     AFTER = 'AAAAAAAAAAAAAAAA'
    #     SEQ = Seq(BEFORE + SITE_SEQ + MID + SITE_SEQ + AFTER, generic_dna)
    #     seq_record = SeqRecord(SEQ)

    #     refactor_context = RefactorContext(seq_record)

    #     occurrences = find_restriction_site_occurrences(
    #             seq_record, RESTRICTION_ENZYME)
    #     self.assertEqual(2, len(occurrences))

    #     seq_record = remove_enzyme_site(
    #             refactor_context, seq_record, RESTRICTION_ENZYME)

    #     occurrences = find_restriction_site_occurrences(
    #             seq_record, RESTRICTION_ENZYME)
    #     self.assertEqual(0, len(occurrences))


    def test_remove_site_in_coding_feature(self):
        """Tests removing a restriction enzyme that falls in a coding
        region.
        """
        RESTRICTION_ENZYME = Restriction.BsmBI
        BEFORE = 'ATGTTTGGGCCCAAATTTGGGAAATTTGGGAAATTTGGGAAATTTGGGAAATTTGGG'
        SITE_SEQ = RESTRICTION_ENZYME.site
        AFTER = 'TAGAAAAAAAAAAAAAAAA'
        SEQ = Seq(BEFORE + SITE_SEQ + AFTER, generic_dna)
        seq_record = SeqRecord(SEQ)

        refactor_context = RefactorContext(seq_record)

        feature_1_loc = FeatureLocation(
                0, len(BEFORE) + len(SITE_SEQ) + 3, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)
        FEATURE_1_SEQ_ORIG = feature_1.extract(str(seq_record.seq))
        FEATURE_1_NUM_CODONS = len(feature_1) / 3

        # Compute fake feature profile.
        fake_profile_values_map = {}
        fake_profile_values_map[feature_1.id] = {
                GCContentFeatureProfile.get_name():
                        [0.2] * FEATURE_1_NUM_CODONS,
                SecondaryStructureFeatureProfile.get_name():
                        [-10] * FEATURE_1_NUM_CODONS,
                CodonRarityFeatureProfile.get_name():
                        [0.5] * FEATURE_1_NUM_CODONS,
        }
        refactor_context.set_feature_id_to_profile_values_map(
                fake_profile_values_map)

        occurrences = find_restriction_site_occurrences(
                seq_record, RESTRICTION_ENZYME)
        self.assertEqual(1, len(occurrences))

        result = _remove_site_in_coding_feature(
                refactor_context, seq_record, occurrences[0], feature_1)

        self.assertTrue(result['is_success'])

        seq_record = result['updated_genome_record']

        FEATURE_1_SEQ_UPDATED = feature_1.extract(str(seq_record.seq))

        occurrences = find_restriction_site_occurrences(
                seq_record, RESTRICTION_ENZYME)
        self.assertEqual(0, len(occurrences))
        self.assertEqual(
                translate_custom(FEATURE_1_SEQ_ORIG),
                translate_custom(FEATURE_1_SEQ_UPDATED))



    def test_remove_site_not_in_coding_feature(self):
        """Test removing a single occurrence of a restriction enzyme.
        """
        RESTRICTION_ENZYME = Restriction.BsmBI
        SITE_SEQ = RESTRICTION_ENZYME.site
        BEFORE = 'CCCCCCCCCCCCCCCCCCCCCCCC'
        AFTER = 'AAAAAAAAAAAAAAAA'
        SEQ = Seq(BEFORE + SITE_SEQ + AFTER, generic_dna)
        seq_record = SeqRecord(SEQ)

        occurrences = find_restriction_site_occurrences(
                seq_record, RESTRICTION_ENZYME)
        self.assertEqual(1, len(occurrences))

        seq_record = _remove_site_not_in_coding_feature(
                seq_record, occurrences[0])

        occurrences = find_restriction_site_occurrences(
                seq_record, RESTRICTION_ENZYME)
        self.assertEqual(0, len(occurrences))


    def test_remove_site_with_partial_feature_overlap__upstream(self):
        """Test fixing a site that is only partially inside a feature.

        We expect a modification to the part that is not inside the feature,
        leaving the feature unchanged.
        """
        RESTRICTION_ENZYME = Restriction.BsmBI
        SITE_SEQ = RESTRICTION_ENZYME.site
        BEFORE = 'CCCCCCCCCCCCCCCCCCCCCCCC'
        AFTER = 'AAAAAAAAAAAAAAAA'
        SEQ = Seq(BEFORE + SITE_SEQ + AFTER, generic_dna)
        seq_record = SeqRecord(SEQ)

        feature_1_loc = FeatureLocation(0, len(BEFORE) + 2, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)
        FEATURE_1_SEQ_ORIG = feature_1.extract(str(seq_record.seq))

        occurrences = find_restriction_site_occurrences(
                seq_record, RESTRICTION_ENZYME)
        self.assertEqual(1, len(occurrences))

        seq_record = _remove_site_with_partial_feature_overlap(
                seq_record, occurrences[0], feature_1)

        FEATURE_1_SEQ_UPDATED = feature_1.extract(str(seq_record.seq))

        occurrences = find_restriction_site_occurrences(
                seq_record, RESTRICTION_ENZYME)
        self.assertEqual(0, len(occurrences))
        self.assertEqual(FEATURE_1_SEQ_ORIG, FEATURE_1_SEQ_UPDATED)


    def test_remove_site_with_partial_feature_overlap__downstream(self):
        """Test fixing a site that is only partially inside a feature.

        We expect a modification to the part that is not inside the feature,
        leaving the feature unchanged.
        """
        RESTRICTION_ENZYME = Restriction.BsmBI
        SITE_SEQ = RESTRICTION_ENZYME.site
        BEFORE = 'CCCCCCCCCCCCCCCCCCCCCCCC'
        AFTER = 'AAAAAAAAAAAAAAAA'
        SEQ = Seq(BEFORE + SITE_SEQ + AFTER, generic_dna)
        seq_record = SeqRecord(SEQ)

        feature_1_loc = FeatureLocation(
                len(BEFORE) + 2, len(SEQ), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)
        FEATURE_1_SEQ_ORIG = feature_1.extract(str(seq_record.seq))

        occurrences = find_restriction_site_occurrences(
                seq_record, RESTRICTION_ENZYME)
        self.assertEqual(1, len(occurrences))

        seq_record = _remove_site_with_partial_feature_overlap(
                seq_record, occurrences[0], feature_1)

        FEATURE_1_SEQ_UPDATED = feature_1.extract(str(seq_record.seq))

        occurrences = find_restriction_site_occurrences(
                seq_record, RESTRICTION_ENZYME)
        self.assertEqual(0, len(occurrences))
        self.assertEqual(FEATURE_1_SEQ_ORIG, FEATURE_1_SEQ_UPDATED)


    # TODO: Figure out why this is failing.
    # def test_find_restriction_occurrences__forward(self):
    #     """Test finding occurences.
    #     """
    #     RESTRICTION_ENZYME = Restriction.BsmBI
    #     SITE_SEQ = RESTRICTION_ENZYME.site
    #     BEFORE = 'CCCCCCCCCCCCCCCCCCCCCCCC'
    #     MID = 'TTTTTTTTTTTT'
    #     AFTER = 'AAAAAAAAAAAAAAAA'
    #     SEQ = Seq(BEFORE + SITE_SEQ + MID + SITE_SEQ + AFTER, generic_dna)
    #     SEQ_RECORD = SeqRecord(SEQ)

    #     OCCURRENCES = find_restriction_site_occurrences(
    #             SEQ_RECORD, RESTRICTION_ENZYME)

    #     # Make sure the occurrences are found.
    #     EXPECTED_INTERVAL_1 = (len(BEFORE), len(BEFORE) + len(SITE_SEQ))
    #     EXPECTED_INTERVAL_2 = (
    #             len(BEFORE) + len(SITE_SEQ) + len(MID),
    #             len(BEFORE) + len(SITE_SEQ) + len(MID) + len(SITE_SEQ))
    #     EXPECTED_OCCURRENCES = [
    #             {
    #                     'interval': EXPECTED_INTERVAL_1,
    #                     'strand': 1
    #             },
    #             {
    #                     'interval': EXPECTED_INTERVAL_2,
    #                     'strand': 1
    #             }
    #     ]
    #     self.assertItemsEqual(EXPECTED_OCCURRENCES, OCCURRENCES)

    #     # Make sure that extracting the sequence gives the actual sequence.
    #     EXTRACTED_SEQ_1 = str(SEQ_RECORD.seq)[
    #             EXPECTED_INTERVAL_1[0]:EXPECTED_INTERVAL_1[1]]
    #     self.assertEqual(SITE_SEQ, EXTRACTED_SEQ_1)
    #     EXTRACTED_SEQ_2 = str(SEQ_RECORD.seq)[
    #             EXPECTED_INTERVAL_2[0]:EXPECTED_INTERVAL_2[1]]
    #     self.assertEqual(SITE_SEQ, EXTRACTED_SEQ_2)


    # TODO: Figure out why this is failing.
    # def test_find_restriction_occurrences__reverse(self):
    #     """Test finding occurences.
    #     """
    #     RESTRICTION_ENZYME = Restriction.BsmBI
    #     SITE_SEQ = reverse_complement(RESTRICTION_ENZYME.site)
    #     BEFORE = 'CCCCCCCCCCCCCCCCCCCCCCCC'
    #     AFTER = 'AAAAAAAAAAAAAAAA'
    #     SEQ = Seq(BEFORE + SITE_SEQ + AFTER, generic_dna)
    #     SEQ_RECORD = SeqRecord(SEQ)

    #     NO_BOUNDS = find_restriction_site_occurrences(
    #             SEQ_RECORD, RESTRICTION_ENZYME)

    #     BOUNDED = find_restriction_site_occurrences(
    #             SEQ_RECORD, RESTRICTION_ENZYME,
    #             start_bound=2, end_bound=len(SEQ) - 2)

    #     equivalent_occurrences = [NO_BOUNDS, BOUNDED]

    #     EXPECTED_INTERVAL = (len(BEFORE), len(BEFORE) + len(SITE_SEQ))
    #     EXPECTED_OCCURRENCES = [{
    #             'interval': EXPECTED_INTERVAL,
    #             'strand': -1
    #     }]
    #     for occurrences in equivalent_occurrences:
    #         # Make sure the occurrences are found.
    #         self.assertItemsEqual(EXPECTED_OCCURRENCES, occurrences)

    #     # Make sure that extracting the sequence gives the actual sequence.
    #     EXTRACTED_SEQ = str(SEQ_RECORD.seq)[
    #             EXPECTED_INTERVAL[0]:EXPECTED_INTERVAL[1]]
    #     self.assertEqual(SITE_SEQ, EXTRACTED_SEQ)


if __name__ == '__main__':
    unittest.main()
