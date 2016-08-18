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
Methods related to calculating and fixing homology issues.
"""

import copy
import difflib
import os

from Bio import SeqIO

from biopython_util import find_features_ending_at
from biopython_util import find_features_starting_at
from biopython_util import get_feature_by_id
from biopython_util import get_genome_record
from biopython_util import InsertType
from biopython_util import update_seq_record_feature
from codon_replacer_util import replace_codons_in_single_feature
from feature_profile import CodonRarityFeatureProfile
from feature_profile import GCContentFeatureProfile
from feature_profile import SecondaryStructureFeatureProfile
from refactor_config import CODONS_TO_REMOVE
from refactor_config import ORIGINAL_CODON_USAGE_MEMEX
from refactor_config import REFACTORED_CODON_USAGE_MEMEX
from refactor_context import RefactorContext

# Where homology matters, we try to get it below this threshold.
HOMOLOGY_THRESHOLD = 0.67

def calc_homology(pair_obj):
    """Calculates the homology between two features:

    The homology is defined as the number of same bases. For simplicity,
    the lengths of the regions are assumed to be the same. So we count up
    the number of matching bases and divide by the length of the sequence.

    Args:
        pair_obj: Object with at least the following keys:
            * source_seq
            * copy_seq

    Returns:
        A float indicating the homology score between the two objects.

    Raises:
        AssertionError if the lengths of the sequences are not the same.
    """
    # Parse the input.
    source_seq = pair_obj['source_seq']
    copy_seq = pair_obj['copy_seq']
    assert len(source_seq) == len(copy_seq)

    # Use the Python difflib module for performing the diff.
    # See: http://docs.python.org/2/library/difflib.html#sequencematcher-objects
    seq_matcher = difflib.SequenceMatcher(a=source_seq, b=copy_seq)
    return seq_matcher.ratio()


def find_features_to_check_for_homology(genome_record):
    """Find the features in the genome_record that we need
    to check for homology.

    Right now, this is a bit hard-coded to our hard overlap fixing strategy.
    Specifically, we look for feature annotations of RBS copies, and then find
    the associated head copy, and then find the corresponding source sequences
    from there.

    TODO: Perhaps we want to create a metadata 'report' while fixing overlaps
        which would keep track of these objects and would be easier to query
        in order to make fixes.

    TODO: This doesn't work for <<<< >>>> type overlaps, but there aren't any
        currently in the region we are dealing with, so figure this out later.

    NOTE: This is wonky for negative strands because of the inconsistent way
        I realized we are annotating the overlapping head copy.

    Returns:
        A list of homology pair objects, with keys:
            * source_id: Feature id of the source.
            * copy_id: Feature id of the copy.
            * source_start_pos: Starting position of the sequence in the source.
            * copy_start_pos: Starting position of the sequence in the copy.
            * source_seq: The actual sequence in the source.
            * copy_seq: The actual sequence in the copy.
    """
    print 'Analyzing copied features for homology...'
    homology_pair_obj_list = []

    rbs_copy_features = filter(
            lambda feature: feature.type == InsertType.FIX_OVERLAP_RBS_COPY,
            genome_record.features
    )

    # From each of these rbs features, we find the associated head copy feature
    # and leverage that to find the copies
    for rbs_copy_feature in rbs_copy_features:
        # print "Reading: " + rbs_copy_feature.id

        rbs_pair_obj = {
                'type': rbs_copy_feature.type,
                'source_id': 'UNKNOWN',
                'source_start': 'UNKNOWN',
                'source_seq': 'UNKNOWN',
                'copy_id': rbs_copy_feature.id,
                'copy_seq': 'UNKNOWN'
        }

        head_copy_pair_obj = {
                'type': 'UNKNOWN',
                'source_id': 'UNKNOWN',
                'source_start': 'UNKNOWN',
                'source_seq': 'UNKNOWN',
                'copy_id': 'UNKNOWN',
                'copy_seq': 'UNKNOWN'
        }

        # Some stuff we can extract without looking at strand, but other stuff
        # we have to do depending on the strand.

        # Get data relevant to the rbs_copy_feature.
        rbs_copy_len = len(rbs_copy_feature)
        rbs_copy_seq = str(rbs_copy_feature.extract(genome_record.seq))
        rbs_pair_obj['copy_seq'] = rbs_copy_seq

        # Get the head feature candidates.
        # NOTE: Might need to do this depending on strand if we fix
        # the way overlapping head copies are annoated. But for now this
        # part where we fetch the candidate works the same for either
        # strand.
        head_copy_feature_candidates = find_features_starting_at(
                rbs_copy_feature.location.end,
                genome_record,
                [InsertType.FIX_OVERLAP_HEAD_COPY]
        )
        assert len(head_copy_feature_candidates) < 2, (
                "Actual len: %d" % len(head_copy_feature_candidates))
        if len(head_copy_feature_candidates) == 1:
            head_copy_feature = head_copy_feature_candidates[0]
            head_copy_len = len(head_copy_feature)
            head_copy_pair_obj['type'] = head_copy_feature.type
            head_copy_pair_obj['copy_id'] = head_copy_feature.id

            # TODO: Dealing with inconsistent overlap head annotation on
            # negative strand. Needs to be cleaned up.
            if rbs_copy_feature.strand == 1:
                head_copy_seq = str(head_copy_feature.extract(genome_record.seq))
            else:
                downstream_feature_candidates = find_features_ending_at(
                    rbs_copy_feature.location.start, genome_record, ['CDS'])
                if len(downstream_feature_candidates) == 1:
                    downstream_feature = downstream_feature_candidates[0]
                    downstream_feature_seq = str(
                            downstream_feature.extract(genome_record.seq))
                    head_copy_seq = str(
                            downstream_feature_seq[:head_copy_len])
                else:
                    continue
            head_copy_pair_obj['copy_seq'] = head_copy_seq

        else:
            # It's unexpected if we don't find the head, break.
            homology_pair_obj_list.append(rbs_pair_obj)
            homology_pair_obj_list.append(head_copy_pair_obj)
            continue

        # TODO: Wonkiness due to the way that we annotated a bit funny. Fix.
        if rbs_copy_feature.strand == 1:
            # Get the source feature these copies are made from.
            source_feature_candidates = find_features_ending_at(
                    rbs_copy_feature.location.start, genome_record, ['CDS'])
        else:
            # Get the source feature these copies are made from.
            source_feature_candidates = find_features_starting_at(
                    rbs_copy_feature.location.end, genome_record, ['CDS'])

        assert len(source_feature_candidates) < 2, (
                "Actual len: %d" % len(source_feature_candidates))
        if len(source_feature_candidates) == 1:
            source_feature = source_feature_candidates[0]
            source_feature_seq = str(source_feature.extract(genome_record.seq))
            rbs_pair_obj['source_id'] = source_feature.id
            rbs_pair_obj['source_start'] = int(source_feature.location.start)
            head_copy_pair_obj['source_id'] = source_feature.id
            head_copy_pair_obj['source_start'] = int(source_feature.location.start)
        else:
            # It's unexpected if we don't find the head, break.
            homology_pair_obj_list.append(rbs_pair_obj)
            homology_pair_obj_list.append(head_copy_pair_obj)
            continue

        # Head source is the last head_copy_len bases of source_feature.
        head_source_seq = source_feature_seq[-1 * head_copy_len:]
        if not head_copy_len == len(head_source_seq):
            homology_pair_obj_list.append(rbs_pair_obj)
            homology_pair_obj_list.append(head_copy_pair_obj)
            continue
        else:
            head_copy_pair_obj['source_seq'] = head_source_seq
            head_copy_pair_obj['source_seq_first_codon_index'] = (
                    (len(source_feature) - head_copy_len) / 3
            )

        # RBS source is the rbs_copy_len bases before head_copy_source.
        rbs_source_seq = source_feature_seq[
                -1 * head_copy_len - rbs_copy_len:
                -1 * head_copy_len]
        if not rbs_copy_len == len(rbs_source_seq):
            homology_pair_obj_list.append(rbs_pair_obj)
            homology_pair_obj_list.append(head_copy_pair_obj)
            continue
        else:
            rbs_pair_obj['source_seq'] = rbs_source_seq
            rbs_pair_obj['source_seq_first_codon_index'] = (
                    (len(source_feature) - head_copy_len - rbs_copy_len) / 3
            )

        # Add the RBS-copy and overlap-head-copy object to the reutrn list.
        homology_pair_obj_list.append(rbs_pair_obj)
        homology_pair_obj_list.append(head_copy_pair_obj)

    return homology_pair_obj_list


def fix_homology_issues(genome_record, ids_to_fix=[]):
    """Finds pairs of copied features created during genome refactoring
    and muddles the upstream original (near the 3' terminus) in order
    to decreate the probability of "snap-back" during insertion.

    Returns:
        A copy of genome_record with homology issues resolved.
    """
    resolved_genome_record = copy.deepcopy(genome_record)

    refactor_context = RefactorContext(resolved_genome_record)

    # Identify features to check for homology issues. These are
    # generally (always?) features that that have head the head/RBS portions
    # copied in order to split apart large overlaps.
    homology_pair_obj_list = find_features_to_check_for_homology(
            resolved_genome_record)

    # Resolve homologies.
    for pair_obj in homology_pair_obj_list:
        copy_id = pair_obj['copy_id']
        if ids_to_fix:
            if not copy_id in ids_to_fix:
                continue
        resolve_single_homology_issue(
                refactor_context,
                pair_obj,
        )

    return resolved_genome_record


def resolve_single_homology_issue(refactor_context, homology_pair_obj):
    """Resolves homology issues between two objects.

    Args:
        refactor_context: The RefactorContext whose SeqRecord will be mutated.
        homology_pair_obj: Object produced by
            find_features_to_check_for_homology(). This is not guaranteed to
            have all the correct keys as of the current implementation.
    """
    genome_record = refactor_context.get_genome_record()

    # The homology_pair_obj passed in may not have all the necessary keys,
    # for example if it's the type of copy we don't know how to deal with yet,
    # so for now we just try-catch any KeyError and report that the homology
    # was not fixed.
    try:
        copy_id = homology_pair_obj['copy_id']
        copy_seq = homology_pair_obj['copy_seq']
        print 'Resolving homology for %s ...' % copy_id

        # Parse the details of the feature to modify from homology_pair_obj.
        feature_to_modify_id = homology_pair_obj['source_id']
        first_codon_to_modify = homology_pair_obj['source_seq_first_codon_index']
        source_feature = get_feature_by_id(genome_record, feature_to_modify_id)
        source_feature_seq = source_feature.extract(genome_record.seq)
        num_codons = len(source_feature) / 3
        avoid_codons_in_positions = {}
        for codon_index in range(first_codon_to_modify, num_codons):
            codon = str(source_feature_seq[codon_index * 3 : codon_index * 3 + 3])
            avoid_codons_in_positions[codon_index] = codon

        # Perform the fix.
        result = replace_codons_in_single_feature(
                refactor_context,
                feature_to_modify_id,
                start_codon_index=first_codon_to_modify,
                avoid_codons_in_positions=avoid_codons_in_positions)
        assert str(result['orig_feature_seq']) != str(result['new_feature_seq'])
        assert result['is_success'], "Resolving homology not successful."

        update_seq_record_feature(
                genome_record,
                feature_to_modify_id,
                result
        )

        # Homology fixed.
        return True

    except KeyError:
        return False

if __name__ == '__main__':
    genome_record = get_genome_record(
            '../data/2013_02_27_22_19_03_mds42_refactored.gbk')
    resolved_genome_record = fix_homology_issues(genome_record)

    # Write the resulting genome.
    genome_output_file = (
            'tmp/temp_fixed.gbk')
    with open(genome_output_file, 'w') as output_fh:
        SeqIO.write(resolved_genome_record, output_fh, 'genbank')
