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
Util methods related to overlap fixing.
"""

from biopython_util import get_region_codon_indeces_in_feature


# Number of bases before a start codon that will capture the RBS.
# When fixing a pair of overlapping genes where the head of a gene is part of
# the overlap, we need to be sure to copy this part of the gene as well.
RBS_BUFFER_SIZE = 20

def does_upstream_feature_overlap_downstream_feature(
        upstream_feature, downstream_feature):
    """Returns a Boolean indicating whether the upstream feature
    overlaps the downstream feature.
    """
    # Answer these questions 2:
    does_upstream_start_before_downstream = (
            upstream_feature.location.start <=
                    downstream_feature.location.start)
    does_downstream_start_before_upstream_ends = (
            downstream_feature.location.start <
                    upstream_feature.location.end)
    return (does_upstream_start_before_downstream and
                    does_downstream_start_before_upstream_ends)


def identify_conflict_region_of_interest(
        upstream_feature, downstream_feature):
    """Returns a tuple pair of genome locations indicating the boundaries
    of the region that we need to check for affected codons. The returned
    interval is "pythonic" in that the first element is inclusive and the
    last element is exclusive.

    For example (ignoring RBS for simplicity):

        xxxxxxx|xxx|
               |xxx|xxxxxx
        0123456|789|

        would return (7,10)

    Region will include:
        * Any overlapping region.
        * RBS region of CDS features.

    If no conflict, returns None.
    """
    # First capture any overlapping region.
    # xxxxxxx|xxx|
    #        |xxx|xxxxxx
    start = downstream_feature.location.start
    end = upstream_feature.location.end

    # Now maybe extend to capture RBS regions.
    # NOTE: Examples use 15 as RBS size, but the constant can be tweaked.
    if downstream_feature.type == 'CDS' and downstream_feature.strand == 1:
        # xxx|xxxxxxxxxxxx|
        #    |         >>>|>>>>>>
        #    |<- 15 ->|
        # Extend the start backward.
        start = downstream_feature.location.start - RBS_BUFFER_SIZE

    if upstream_feature.type == 'CDS' and upstream_feature.strand == -1:
        #              |<- 15 ->|
        # <<<<<<<<<|<<<         |
        #          |xxxxxxxxxxxx|xxxxxx
        # Extend the end backward.
        end = upstream_feature.location.end + RBS_BUFFER_SIZE

    # Check that the region found is none-trivial and return.
    if end > start:
        return (int(start), int(end))

    # No valid conflicting region found.
    return None


def is_feature_in_region(feature, region):
    """Calculate whether the feature overlaps the region.
    """
    return (region[0] < feature.location.end and
            region[1] > feature.location.start)


def does_feature_have_forbidden_codon_in_region(
        feature, region, seq_record, forbidden_codons):
    """Checks the feature codons that fall in the region.
    """
    return does_feature_have_codon_in_list_region(
            feature, region, seq_record, forbidden_codons)


def does_feature_have_codon_in_list_region(
        feature, region, seq_record, codon_list, return_codon_index=False):
    """Checks the feature codons that fall in the region against
    the codon_list. If any match, returns True. Else False.
    """
    feature_seq = str(feature.extract(seq_record.seq))

    codon_indeces = get_region_codon_indeces_in_feature(
            feature, region)

    # Iterate through the codons, checking for any that are in target list.
    for codon_index in codon_indeces:
        codon = feature_seq[codon_index * 3 : codon_index * 3 + 3]
        if codon in codon_list:
            if return_codon_index:
                return codon_index
            else:
                return True

    # No forbidden codons found.
    return False
