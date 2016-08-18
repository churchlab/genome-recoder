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
Tests for gc_content_fixer.py
"""

import copy
import random
import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from gc_content_fixer import automated_intergenic_gc_fixer
from gc_content_fixer import _calculate_feature_annotation_shadow
from gc_content_fixer import GCContentConstraints


class TestAutomatedIntergenicGcFixer(unittest.TestCase):

    def test_noop__no_intervals(self):
        seq = Seq(''.join([random.choice('ATGC') for i in range(200)]),
                generic_dna)
        orig_seq_record = SeqRecord(seq)
        seq_record = copy.deepcopy(orig_seq_record)
        interval_list = []
        automated_intergenic_gc_fixer(seq_record, interval_list)
        self.assertEqual(str(orig_seq_record.seq), str(seq_record.seq))

    def test_noop__small_interval(self):
        # Sequence is all GC so no changes to make.
        seq = Seq(''.join([random.choice('GC') for i in range(200)]),
                generic_dna)
        orig_seq_record = SeqRecord(seq)
        seq_record = copy.deepcopy(orig_seq_record)
        interval_list = [(100, 125)]
        automated_intergenic_gc_fixer(seq_record, interval_list)
        self.assertEqual(str(orig_seq_record.seq), str(seq_record.seq))

    def test_change_some(self):
        """Keep changing positions until target is met.
        """
        # Repeat to make sure not stochastic result
        # ...though it still could be :)
        for i in range(100):
            # Sequence is all AT so some will have to change.
            seq = Seq(''.join([random.choice('AT') for i in range(200)]),
                    generic_dna)
            orig_seq_record = SeqRecord(seq)
            seq_record = copy.deepcopy(orig_seq_record)
            self.assertTrue(GC(seq_record.seq) == 0)
            interval_list = [(100, 125)]
            automated_intergenic_gc_fixer(seq_record, interval_list)
            for center in range(100, 125):
                window_seq = seq_record.seq[center - 50:center + 50]
                self.assertTrue(GC(window_seq) >= 30,
                        "GC for window centered at %d is %f" % (
                                center, GC(window_seq)))
            self.assertNotEqual(str(orig_seq_record.seq), str(seq_record.seq))

    def test_change_all(self):
        """Change all positions around intervals.
        """
        # Sequence is all AT so some will have to change.
        seq = Seq(''.join([random.choice('AT') for i in range(200)]),
                generic_dna)
        orig_seq_record = SeqRecord(seq)
        seq_record = copy.deepcopy(orig_seq_record)
        self.assertTrue(GC(seq_record.seq) < 30)
        interval_list = [(100, 125)]
        constraint_obj = GCContentConstraints()
        constraint_obj.local_window_lower_bound = 1.0
        constraint_obj.local_window_upper_bound = 1.1 # no upper bound
        automated_intergenic_gc_fixer(seq_record, interval_list,
                gc_content_constraint_obj=constraint_obj)
        for center in range(100, 125):
            window_seq = seq_record.seq[center - 50:center + 50]
            self.assertTrue(GC(window_seq) == 100,
                    "GC for window centered at %d is %f" % (
                            center, GC(window_seq)))
        self.assertNotEqual(str(orig_seq_record.seq), str(seq_record.seq))

    def test_change_all__multiple_intervals(self):
        """Change all positions around intervals.
        """
        # Sequence is all AT so some will have to change.
        seq = Seq(''.join([random.choice('AT') for i in range(300)]),
                generic_dna)
        orig_seq_record = SeqRecord(seq)
        seq_record = copy.deepcopy(orig_seq_record)
        self.assertTrue(GC(seq_record.seq) < 30)
        interval_list = [(100, 125), (160, 180)]
        constraint_obj = GCContentConstraints()
        constraint_obj.local_window_lower_bound = 1.0
        constraint_obj.local_window_upper_bound = 1.1 # no upper bound
        automated_intergenic_gc_fixer(seq_record, interval_list,
                gc_content_constraint_obj=constraint_obj)
        for interval in interval_list:
            for center in range(*interval):
                window_seq = seq_record.seq[center - 50:center + 50]
                self.assertTrue(GC(window_seq) == 100,
                        "GC for window centered at %d is %f" % (
                                center, GC(window_seq)))
        self.assertNotEqual(str(orig_seq_record.seq), str(seq_record.seq))

    def test_calculate_feature_annotation_shadow(self):
        # Sequence is all AT so some will have to change.
        seq = Seq(''.join([random.choice('AT') for i in range(200)]),
                generic_dna)
        seq_record = SeqRecord(seq)

        # We will expect a 60 base shadow from 60 - 120. This is the CDS
        # and 20 bases upstream.
        feature_1_loc = FeatureLocation(80, 120, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        # Calculate the shadow.
        pos_blacklist_set = _calculate_feature_annotation_shadow(seq_record)

        # Check it.
        EXPECTED_BLACKLIST_SET = set(list(range(60, 120)))
        self.assertEqual(EXPECTED_BLACKLIST_SET, pos_blacklist_set)

    def test_avoid_changes_in_shadows(self):
        """Avoid changing bases in CDS features.
        """
        # Sequence is all AT so some will have to change.
        seq = Seq(''.join([random.choice('AT') for i in range(200)]),
                generic_dna)
        orig_seq_record = SeqRecord(seq)

        # We will expect a 60 base shadow from 60 - 120. This is the CDS
        # and 20 bases upstream.
        feature_1_loc = FeatureLocation(80, 120, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        orig_seq_record.features.append(feature_1)

        seq_record = copy.deepcopy(orig_seq_record)

        self.assertTrue(GC(seq_record.seq) == 0)

        # Hit just one position.
        interval_list = [(100, 101)]

        # Aim for 100% to hit all bases possible.
        constraint_obj = GCContentConstraints()
        constraint_obj.local_window_lower_bound = 1.0
        constraint_obj.local_window_upper_bound = 1.1 # no upper bound
        automated_intergenic_gc_fixer(seq_record, interval_list,
                gc_content_constraint_obj=constraint_obj)

        # Expect window centered at 100 to have GC 60, which is all
        window_seq = seq_record.seq[50:150]
        self.assertEqual(40, GC(window_seq))

        # Make sure shadow seq is unchanged.
        self.assertEqual(
            str(feature_1.extract(orig_seq_record.seq)),
            str(feature_1.extract(seq_record.seq)))


if __name__ == '__main__':
    unittest.main()
