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
Tests for codon_replacer_util.py.
"""

import os
import unittest

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from codon_replacer_util import get_essential_feature_ids
from codon_replacer_util import replace_codons_in_feature_subset
from codon_replacer_util import replace_codons_in_single_feature
from codon_replacer_util import replace_forbidden_codons
from codon_replacer_util import update_seq_record_with_partial_results
from feature_profile import GCContentFeatureProfile
from feature_profile import SecondaryStructureFeatureProfile
from feature_profile import CodonRarityFeatureProfile
from refactor_config import ORIGINAL_CODON_USAGE_MEMEX
from refactor_config import REFACTORED_CODON_USAGE_MEMEX
from refactor_config import ORIGINAL_GENOME_RECORD
from refactor_context import RefactorContext


class TestUtilMethods(unittest.TestCase):
    """Test exported utility methods.

    Most of these simply test calling the methods. We should add tests that
    actually try the logic as time allows.
    """

    def test_get_essential_feature_ids(self):
        get_essential_feature_ids(SeqRecord(Seq('ATG')))


    def test_replace_codons_in_feature_subset(self):
        TMP_RESULT_FILE = 'tmp_test_replace_codons_in_feature_subset'
        replace_codons_in_feature_subset(
                SeqRecord(Seq('ATG')),
                [], # essential_feature_ids
                [], # codons_to_remove
                ORIGINAL_CODON_USAGE_MEMEX,
                REFACTORED_CODON_USAGE_MEMEX,
                {}, # feature_id_to_profile_values_map
                0, # range_start
                0, # range_end
                TMP_RESULT_FILE,
                False) # debug

        if os.path.exists(TMP_RESULT_FILE):
            os.remove(TMP_RESULT_FILE)


    def test_replace_codons_in_single_feature(self):
        FEATURE_ID = 'ECMDS42_0001'
        refactor_context = RefactorContext(ORIGINAL_GENOME_RECORD)
        refactor_context.set_feature_id_to_profile_values_map({
            FEATURE_ID: {
                GCContentFeatureProfile.get_name(): [0.375, 0.325, 0.35, 0.325, 0.35, 0.4, 0.45, 0.425, 0.45, 0.475, 0.475, 0.475, 0.5, 0.475, 0.5, 0.475, 0.5, 0.5, 0.5, 0.525, 0.575, 0.575],

                SecondaryStructureFeatureProfile.get_name(): [-0.6, -0.6, -0.6, -1.3, -1.3, -2.5, -2.5, -2.5, -2.3, -1.9, -1.9, -1.9, -1.9, -1.9, -1.9, 2.1, -1.5, -3.7, -7.4, -9.4, -12.0, -12.5],

                CodonRarityFeatureProfile.get_name(): [0.368, 0.434, 0.488, 0.5640000000000001, 0.5900000000000001, 0.47800000000000004, 0.414, 0.43599999999999994, 0.422, 0.454, 0.454, 0.45, 0.43600000000000005, 0.45, 0.45, 0.388, 0.372, 0.394, 0.36000000000000004, 0.34400000000000003, 0.4000000000000001, 0.39]

        }})

        replace_codons_in_single_feature(refactor_context, FEATURE_ID)


    def test_replace_forbidden_codons(self):
        GENOME_RECORD = SeqRecord(Seq('ATG'))
        refactor_context = RefactorContext(GENOME_RECORD)
        refactor_context.set_feature_id_to_profile_values_map({})
        replace_forbidden_codons(refactor_context)


    def test_replace__debugging_forbiddenCodon_annotation(self):
        ECMDS42_3701_SEQ = 'ATGCAACCTTTTGGCGTACTTGACCGCTATATCGGTAAAACTATTTTCACCACCATCATGATGACACTGTTCATGCTGGTGTCGCTGTCGGGCATTATCAAGTTTGTCGATCAGCTGAAAAAAGCCGGGCAGGGGAGTTACGACGCGTTAGGCGCAGGAATGTATACCTTGCTGAGCGTGCCGAAAGATGTGCAGATCTTCTTCCCGATGGCGGCTCTGCTTGGGGCGTTGCTTGGTCTTGGGATGCTGGCGCAGCGCAGCGAACTGGTGGTGATGCAGGCTTCTGGTTTTACCCGTATGCAGGTGGCGCTGTCGGTGATGAAAACCGCCATTCCGCTGGTCTTGCTGACGATGGCGATTGGCGAATGGGTCGCGCCGCAGGGCGAGCAGATGGCGCGTAACTACCGTGCGCAGGCGATGTACGGCGGCTCGTTGCTCTCTACCCAGCAAGGCTTATGGGCGAAAGATGGCAACAACTTCGTCTACATTGAGCGGGTTAAAGGTGACGAAGAGTTAGGTGGCATCAGCATTTATGCCTTTAACGAGAATCGTCGTCTGCAATCCGTACGCTATGCCGCTACTGCGAAGTTTGACCCGGAACATAAAGTCTGGCGTCTGTCGCAGGTTGATGAATCTGATCTGACCAATCCGAAACAGATTACCGGTTCGCAGACGGTGAGCGGCACCTGGAAAACCAACCTCACGCCGGACAAACTGGGCGTGGTGGCGCTGGACCCGGATGCACTCTCTATCAGCGGTTTGCACAACTATGTGAAGTATCTGAAGTCGAGCGGTCAGGATGCCGGACGTTATCAGCTCAACATGTGGAGCAAAATCTTCCAGCCGCTATCTGTGGCGGTGATGATGCTGATGGCGCTGTCGTTCATCTTTGGCCCACTGCGTAGCGTACCGATGGGCGTGCGTGTGGTCACCGGTATCAGTTTCGGTTTTGTCTTCTACGTACTGGACCAGATCTTCGGCCCGCTGACGTTGGTTTATGGCATCCCGCCGATCATCGGCGCACTGTTGCCAAGCGCCAGCTTCTTCTTAATCAGCCTGTGGCTGTTAATGAGAAAATCGTAA'
        seq_record = SeqRecord(ECMDS42_3701_SEQ)

        feature_1_loc = FeatureLocation(0, len(ECMDS42_3701_SEQ), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='ECMDS42_3701')
        seq_record.features.append(feature_1)

        REFACTOR_CONTEXT = RefactorContext(seq_record)

        replace_codons_in_single_feature(REFACTOR_CONTEXT, feature_1.id)


if __name__ == '__main__':
    unittest.main()
