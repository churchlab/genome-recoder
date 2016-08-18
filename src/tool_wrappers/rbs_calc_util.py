"""
Module that adapts the Salis rbs_calc package to the way that we use it.
"""

import math
import os
import pickle

from Bio.Seq import reverse_complement

from biopython_util import get_genome_record
from rbs_calc.RBS_Calculator import RBS_Calculator


CACHE_DIR = 'cache'

CACHE_LOCATION = os.path.join(CACHE_DIR, 'mds42_rbs_profile.cache')

MDS42_SOURCE_GENBANK = '../data/mds42_full.gbk'

FEATURES_TO_PROFILE_SET = set(['CDS', 'misc_feature'])

RBS_CALC_BUFFER = 30

UNKNOWN_RBS_CALC_RESULT = 'UNKNOWN'


def get_mds42_rbs_strength_profile(force_recalculate=False):
    """Returns the RBS strength profile for each feature in the genome,
    using cached results, or calculating again if necessary.

    Returns:
        Dictionary mapping locus_tag to RBS strength.
    """
    print 'Getting MDS42 RBS strength profile...'

    if not force_recalculate:
        # First try to get the results from cache.
        try:
            with open(CACHE_LOCATION) as fh:
                print '...Using cached RBS profile from %s' % CACHE_LOCATION
                return pickle.load(fh)
        except:
            # Most likely cache doesn't exist.
            pass

    genome_record = get_genome_record(MDS42_SOURCE_GENBANK)
    forward_seq = genome_record.seq

    features_to_profile = [feature for feature in genome_record.features
            if feature.type in FEATURES_TO_PROFILE_SET]

    locus_tag_to_rbs_strength_map = {}

    num_features = len(features_to_profile)
    for feature_index in range(num_features):
        feature = features_to_profile[feature_index]

        print 'Calculating feature %d of num_features %d' % (
                    feature_index + 1, num_features)

        result = calc_rbs_score_for_feature(
                feature, forward_seq, RBS_CALC_BUFFER)
        locus_tag = result['locus_tag']
        expression = result['expression']
        locus_tag_to_rbs_strength_map[locus_tag] = expression

    # Cache the results.
    if not os.path.exists(CACHE_DIR):
        os.mkdir(CACHE_DIR)
    with open(CACHE_LOCATION, 'w') as fh:
        pickle.dump(locus_tag_to_rbs_strength_map, fh)

    return locus_tag_to_rbs_strength_map


def calc_rbs_score_for_feature(feature, forward_seq,
        rbs_calc_buffer=RBS_CALC_BUFFER):
    """Calculate the rbs score for a single feature.
    """
    # Figure out the polarity-aware seq.
    if feature.strand == 1:
        seq = str(forward_seq[
                feature.location.start - rbs_calc_buffer:
                feature.location.start + rbs_calc_buffer])
    else:
        forward_read = forward_seq[
                feature.location.end - rbs_calc_buffer:
                feature.location.end + rbs_calc_buffer]
        seq = str(reverse_complement(forward_read))

    assert 2 * rbs_calc_buffer == len(seq), (
            "Unexpected sequence length for RBS calc. Debug.")

    # Perform calculation.
    start_range = (rbs_calc_buffer, rbs_calc_buffer)
    return do_calc_rbs(feature, seq, start_range, rbs_calc_buffer)


def calc_internal_rbs_expression(seq, accept_rbs_after_pos=0):
    """Attempts to calculate internal RBS expression for the given sequence.
    """
    start_range = [accept_rbs_after_pos, len(seq)]
    rbs_calc_obj = RBS_Calculator(seq, start_range)
    rbs_calc_obj.calc_dG()

    # List corresponding to the length of the number of start codons found
    # in the mRNA of interest.
    dG_total_list = rbs_calc_obj.dG_total_list[:]

    # Finish parsing the rsults.
    if dG_total_list:
        start_pos_list = rbs_calc_obj.start_pos_list[:]
        kinetic_score_list = rbs_calc_obj.kinetic_score_list[:]
        putative_start = start_pos_list[0]
        dG = dG_total_list[0]
        kinetic_score = kinetic_score_list[0]
        expression = rbs_calc_obj.calc_expression_level(dG)
    else:
        putative_start = UNKNOWN_RBS_CALC_RESULT
        dG = UNKNOWN_RBS_CALC_RESULT
        kinetic_score = UNKNOWN_RBS_CALC_RESULT
        expression = UNKNOWN_RBS_CALC_RESULT

    result = {
        'dG': dG,
        'kinetic_score': kinetic_score,
        'expression': expression,
        'seq': seq,
        'putative_start': putative_start
    }
    return result


def do_calc_rbs(feature, seq, start_range, rbs_calc_buffer=None):
    """Light wrapper helper method that makes the call and parses the results
    given the exact seqence to send.
    """
    rbs_calc_obj = RBS_Calculator(seq, start_range)
    rbs_calc_obj.calc_dG()

    # List corresponding to the length of the number of start codons found
    # in the mRNA of interest.
    dG_total_list = rbs_calc_obj.dG_total_list[:]

    # The above fails in MDS42 for some misc_features, but for CDS only when
    # the start codon is ATT, which Salis's calculator doesn't handle.
    # So attempt to remedy this case by temporarily swapping in an ATG
    # at that position.
    if (not rbs_calc_buffer is None and not dG_total_list and
            ('ATT' == seq[rbs_calc_buffer : rbs_calc_buffer + 3])):
        mod_seq = (
                seq[:rbs_calc_buffer] + 'ATG' + seq[rbs_calc_buffer + 3:]
        )
        assert len(seq) == len(mod_seq)

        rbs_calc_obj = RBS_Calculator(mod_seq, start_range)
        rbs_calc_obj.calc_dG()
        dG_total_list = rbs_calc_obj.dG_total_list[:]

    # Finish parsing the rsults.
    if dG_total_list:
        start_pos_list = rbs_calc_obj.start_pos_list[:]
        kinetic_score_list = rbs_calc_obj.kinetic_score_list[:]
        if not rbs_calc_buffer is None:
            assert start_pos_list[0] == rbs_calc_buffer
            assert len(start_pos_list) == 1
        dG = dG_total_list[0]
        kinetic_score = kinetic_score_list[0]
        expression = rbs_calc_obj.calc_expression_level(dG)
    else:
        dG = UNKNOWN_RBS_CALC_RESULT
        kinetic_score = UNKNOWN_RBS_CALC_RESULT
        expression = UNKNOWN_RBS_CALC_RESULT

    result = {
        'dG': dG,
        'kinetic_score': kinetic_score,
        'expression': expression,
        'seq': seq
    }

    # HACK(gleb, 6/4/14): Allow feature=None for analyzing AGR.
    if not feature is None:
        result.update({
            'locus_tag':
                feature.qualifiers.get('locus_tag', ['UNKNOWN_LOCUS'])[0],
            'gene':
                feature.qualifiers.get('gene', ['UNKNOWN_GENE'])[0],
        })


    return result


def is_rbs_strength_conserved(
        orig_expression_profile,
        feature,
        underlying_seq,
        allowed_order_of_magnitude_error):
    """Checks whether the RBS expression level as measured by the Salis
    tool is within a couple orders of magnitude of the features original
    expression.

    Args:
        orig_expression_profile: Original mapping from locus_tag to expression.
        feature: The feature we're checking. Most relevant are the feature
            id and location.
        underlying_seq: The underlying sequence that the locations will
            be applied against. So clients must figure this out and pass
            it depending on the context.
    Returns:
        A Boolean indicating whether RBS is conserved.
    """
    # Get the original expression from cache.
    orig_expression = orig_expression_profile[feature.id]

    # Calculate the new expression.
    new_expression = calc_rbs_score_for_feature(
            feature, underlying_seq, RBS_CALC_BUFFER)['expression']
    if new_expression == UNKNOWN_RBS_CALC_RESULT:
        # The only case this is known to happen is if the start codon
        # is ATT, so we temporariliy swap out ATG for the calc.
        orig_feature_seq = str(underlying_seq.extract(feature))
        assert 'ATT' == orig_feature_seq[0:3]

        # Modify the underlying sequence.
        if feature.strand == 1:
            mod_seq = (
                    underlying_seq[0:feature.location.start] +
                    'ATG' +
                    underlying_seq[feature.location.start:]
            )
        else:
            mod_seq = (
                    underlying_seq[:feature.location.end - 3] +
                    'CAT' +
                    underlying_seq[feature.location.end:]
            )
        assert len(underlying_seq) == len(mod_seq)

        # Recalculate expression.
        new_expression = calc_rbs_score_for_feature(
                feature, mod_seq, RBS_CALC_BUFFER)['expression']
        assert not isinstance(new_expression, str)

    # Determine whether the log expression is within the allowed error order
    # of magnitude.
    log_orig_expression = math.log10(orig_expression)
    log_new_expression = math.log10(new_expression)
    return (abs(log_new_expression - log_orig_expression) <
            allowed_order_of_magnitude_error)


if __name__ == '__main__':
    get_mds42_rbs_strength_profile(force_recalculate=True)
