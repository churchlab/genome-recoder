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
Methods for checking that Gen9 synthesis constraints are satisfied.
"""

import re

from Bio.SeqUtils import GC

from homopolymers import HOMOPOLYMER_RE
from refactor_config import RESTRICTION_ENZYME_SITES_TO_REMOVE
from restriction_sites import find_all_restriction_sites_in_seq


def find_homopolymer_runs(seq):
    return [(i.start(), i.group()) for i in re.finditer(HOMOPOLYMER_RE, seq)]


def check_all(seq):
    return {
        'homopolymers': find_homopolymer_runs(seq),
        'restriction_sites': find_all_restriction_sites_in_seq(seq,
                RESTRICTION_ENZYME_SITES_TO_REMOVE)
    }


def analyze_seq_GC(seq, window_size=100, gc_lower=30, gc_upper=75):
    """Rudimentary method that prints regions of bad GC content.
    """
    half_window_size = window_size / 2

    # Identify bounds for which we can center a window seq (half the
    # window size from each end).
    seq_lower_bound = half_window_size
    seq_upper_bound = len(seq) - half_window_size

    # Store data for this sub_fragment in this list here.
    for window_center in range(len(seq)):
        if seq_lower_bound < window_center < seq_upper_bound:
            window_start = window_center - half_window_size
            window_end = window_center + half_window_size
            window_seq = seq[window_start:window_end]
            assert len(window_seq) == window_size, (
                    'Actual: %d' % len(window_seq))
            gc = GC(window_seq)
            if gc < gc_lower or gc > gc_upper:
                print 'VIOLATION', gc, window_start, window_seq
