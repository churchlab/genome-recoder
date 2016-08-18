"""
Methods for finding pairs of conflicting genes.

NOTE: This is an effort to unclutter overlaps.py where the logic for both
finding and fixing overlaps is currently located.
"""

import os
import pickle

from conflicting_pair_common import does_feature_have_forbidden_codon_in_region
from conflicting_pair_common import does_upstream_feature_overlap_downstream_feature
from conflicting_pair_common import identify_conflict_region_of_interest
from conflicting_pair_common import is_feature_in_region
from conflicting_pair_common import RBS_BUFFER_SIZE
from refactor_config import CACHE_DIR


FIND_OVERLAPS_CACHE_LOCATION = os.path.join(
        CACHE_DIR, 'mds42_find_overlaps.cache')


def find_all_overlaps(
        genome_record,
        forbidden_codons,
        require_feature_type=[],
        ignore_feature_type=[],
        essential_feature_ids=[],
        full_analysis=False,
        cache=False,
        include_close_features=False):
    """Finds overlaps among all features of a Genome.

    Args:
        genome_record: A SeqRecord representing the genome
            and all of its features, as read, say, from a genbank
            file.
        require_feature_type: List of feature types to limit overlap
            checking to (e.g. 'CDS').
        ignore_feature_type: List of feature types to ignore while checking
            overlaps (e.g. 'gene', since it seems there is a 'gene' for
            every 'CDS').
        essential_feature_ids: List of features that we care about. If this
            is provided, the require and ignore types above are ignored.
        full_analysis: Whether to compute metadata for overlaps. Useful for
            analysis, but not necessarily when getting refactoring done.
        cache: Whether to use and store cached results.  This is useful in
            when debugging fixing overlaps, and not having to find all overlaps
            over and over.
        include_close_features: If True, broaden the check to look for features
            that are close enough to potentially interfere with either or both
            RBS regions.

    Returns:
        A list of dictionary objects with keys:
            * upstream_feature
            * downstream_feature

            and if full_analysis=True, additional metadata keys:
                * overlap_type
                * overlap_size
                * overlap_polarity
    """
    print 'Finding all overlaps ...'

    # First try to recover cached results:
    if cache:
        try:
            with open(FIND_OVERLAPS_CACHE_LOCATION) as fh:
                overlaps_found = pickle.load(fh)
                print 'Using cached results from %s' % (
                        FIND_OVERLAPS_CACHE_LOCATION,)
                return overlaps_found
        except:
            # Most likely cache doesn't exist.
            pass


    # Without loss of generality, we find overlaps where the first feature is
    # upstream of the other feature.

    overlaps_found = []

    ### Identify the features we are going to check.
    features_to_check = genome_record.features

    # If essential_feature_ids is present, the others are ignored.
    if essential_feature_ids:
        essential_feature_id_set = set(essential_feature_ids)
        features_to_check = filter(
                lambda feature: feature.id in essential_feature_id_set,
                features_to_check)
    else:
        if len(require_feature_type) > 0:
            features_to_check = filter(
                    lambda feature: feature.type in require_feature_type,
                    features_to_check)
        if len(ignore_feature_type) > 0:
             features_to_check = filter(
                    lambda feature: not feature.type in ignore_feature_type,
                    features_to_check)

    # First, order the features according to start position so that we can do
    # less then n^2 comparisons.
    sorted_features = sorted(features_to_check,
            cmp=lambda x,y: cmp(x.location.start, y.location.start))

    num_features = len(sorted_features)
    for upstream_feature_index in range(0, num_features):
        for downstream_feature_index in range(
                upstream_feature_index + 1, num_features):
            upstream_feature = sorted_features[upstream_feature_index]
            downstream_feature = sorted_features[downstream_feature_index]

            # Perform the desired overlap check (include_close_features)
            # is the more recent strategy, breaking out of loop if not found,
            # otherwise populating the conflict variable.
            if include_close_features:
                conflict_check = (
                        is_upstream_downstream_pair_candidate_for_separation(
                                upstream_feature, downstream_feature))

                if conflict_check['is_conflict']:
                    conflict = {
                            'upstream_feature': upstream_feature,
                            'downstream_feature': downstream_feature,
                            'conflict_type': conflict_check['conflict_type']
                    }
                else:
                    continue
            else:
                is_overlap = does_upstream_feature_overlap_downstream_feature(
                        upstream_feature, downstream_feature)
                if is_overlap:
                    conflict = {
                            'upstream_feature': upstream_feature,
                            'downstream_feature': downstream_feature,
                            'conflict_type': 'overlap'
                    }
                else:
                    continue

            # Now do a second level of checks that see whether the conflict
            # needs to be fixed. For example, if there are no synonymous codons
            # in the afflicted regions, it's not necessary to fix the conflict.
            need_to_fix_conflict = does_conflict_need_fixing(
                    conflict, genome_record, forbidden_codons)

            if not need_to_fix_conflict:
                continue

            assert conflict, "No conflict, but came here. That's a bug."

            # Special overlap type handling.
            if conflict['conflict_type'] == 'overlap':
                # Determine type and size.
                if (downstream_feature.location.end <=
                        upstream_feature.location.end):
                    overlap_type = 'upstream_contains_downstream'
                    overlap_size = (downstream_feature.location.end -
                        downstream_feature.location.start)
                else:
                    overlap_type = 'misc'
                    overlap_size = (upstream_feature.location.end -
                        downstream_feature.location.start)

                conflict['overlap_type'] = overlap_type
                conflict['overlap_size'] = overlap_size

            if full_analysis:
                # Determine polarity.
                overlap_polarity = get_overlap_polarity_string(
                        upstream_feature, downstream_feature)
                conflict['polarity'] = overlap_polarity

            # Add to overlaps list and continue looping.
            overlaps_found.append(conflict)

    # Maybe cache results.
    if cache:
        if not os.path.exists(CACHE_DIR):
            os.mkdir(CACHE_DIR)
        with open(FIND_OVERLAPS_CACHE_LOCATION, 'w') as fh:
            pickle.dump(overlaps_found, fh)

    print 'Found ' + str(len(overlaps_found)) + ' overlaps'
    return overlaps_found


def is_upstream_downstream_pair_candidate_for_separation(
        upstream_feature, downstream_feature):
    """Determines whether a pair of features is a candidate for
    separation analysis. This is the case if the features overlap,
    or if one of the features potentially has its RBS region within
    the other feature.

    The upstream feature must start at or before the downstream feature.

    Returns:
        An object with keys:
            * is_conflict: {Boolean}
            * conflict_type: {String}
    """
    is_upstream_start_same_or_before_downstream = (
            upstream_feature.location.start <=
                    downstream_feature.location.start)
    assert is_upstream_start_same_or_before_downstream

    # Obvious case is if they overlap.
    is_overlap = does_upstream_feature_overlap_downstream_feature(
            upstream_feature, downstream_feature)
    if is_overlap:
        return {
            'is_conflict': True,
            'conflict_type': 'overlap'
        }

    # Less obvious is if they are too close and may affect each other's
    # RBS.
    maybe_rbs_conflict = is_potential_rbs_conflict(
            upstream_feature, downstream_feature)
    if maybe_rbs_conflict:
        return {
            'is_conflict': True,
            'conflict_type': 'rbs_conflict'
        }

    # All checks passed, these features are not candidates for separation.
    return {
        'is_conflict': False,
        'conflict_type': None
    }


def is_potential_rbs_conflict(
        upstream_feature, downstream_feature):
    """Check if the features are close enough to affect either or both
    RBS features.

    If the 5' end of either of the feature is within RBS_BUFFER_SIZE
    bp's from the other feature, unless they are
    facing each other in which case they would not affect each other's RBS.
    """
    # NOTE: This could be written as fancy boolean combos, but a gauntlet of
    # if is more clear to understand.

    # Both must be CDS, since RBS for one will only be affected if
    # we accidentally break it while doing synonymous changes to the other,
    # but we only do synonymous swaps in CDS regions.
    are_both_cds = (
            upstream_feature.type == 'CDS' and
            downstream_feature.type == 'CDS'
    )
    if not are_both_cds:
        return False

    # Check that they are close enough.
    is_distance_small_enough = (
            downstream_feature.location.start - upstream_feature.location.end <
                    RBS_BUFFER_SIZE)
    if not is_distance_small_enough:
        return False

    # Case where they are head on doesn't matter.
    # >>>>>>>>>> <<<<<<<<<<
    are_head_on = (upstream_feature.strand == 1 and
            downstream_feature.strand == -1)
    if are_head_on:
        return False

    # Made it through the gauntlet so this is a potential rbs conflict.
    return True


def does_conflict_need_fixing(conflict_obj, seq_record, forbidden_codons):
    """
    Args:
        conflict_obj: Object with keys:
            * upstream_feature
            * downstream_feature
            * conflict_type
        seq_record: The genome record containing the underlying sequence.
        forbidden_codons: Set of codons that we are removing.
    Returns:
        Boolean indicating whether the conflict needs to be fixed or not.
    """
    # Parse the input.
    upstream_feature = conflict_obj['upstream_feature']
    downstream_feature = conflict_obj['downstream_feature']
    conflict_type = conflict_obj['conflict_type']

    # If not an overlap, then both must be CDS in order for there to even be a
    # possible conflict, as issues with none-overlaps only arise if one or
    # the other CDS's RBS is mangled during synonymous codon replacement.
    if (not conflict_type == 'overlap' and
            not (upstream_feature.type == 'CDS' and
                    downstream_feature.type == 'CDS')):
        return False

    # Now, identify the relevant region of the genome to look at. That is,
    # we'll look in the affected features in this region.
    region_of_interest = identify_conflict_region_of_interest(
            upstream_feature, downstream_feature)
    if not region_of_interest:
        return False

    # Check whether either feature is a CDS with forbidden codons that
    # intersect the region of interest.
    # This is a fairly cheap check so we can do it for any size region.
    upstream_has_forbidden = (
            upstream_feature.type == 'CDS' and
            is_feature_in_region(upstream_feature, region_of_interest) and
            does_feature_have_forbidden_codon_in_region(
                    upstream_feature, region_of_interest, seq_record,
                    forbidden_codons)
    )
    downstream_has_forbidden = (
            downstream_feature.type == 'CDS' and
            is_feature_in_region(downstream_feature, region_of_interest) and
            does_feature_have_forbidden_codon_in_region(
                    downstream_feature, region_of_interest, seq_record,
                    forbidden_codons)
    )
    if not upstream_has_forbidden and not downstream_has_forbidden:
        # No forbidden codons in region, nothing to fix.
        return False

    # Since we got here, it means there is a potential issue to fix.
    return True


def get_overlap_polarity_string(feature_1, feature_2):
    if feature_1.strand == 1 and feature_2.strand == 1:
        return '>>'
    elif feature_1.strand == 1 and feature_2.strand == -1:
        return '><'
    elif feature_1.strand == -1 and feature_2.strand == 1:
        return '<>'
    else:
        return '<<'
