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
Methods for checking the results of refactoring for mistakes.
"""

from Bio.Data.CodonTable import TranslationError

from biopython_util import does_feature_start_with_valid_start_codon
from biopython_util import does_interval_overlap_feature
from biopython_util import get_feature_by_id
from biopython_util import translate_custom
from refactor_config import CODONS_TO_REMOVE
from refactor_config import SELENOCYSTEINE_CODON
from refactor_config import TRANSLATION_CONSERVED_EXCEPTIONS
from tool_wrappers.rbs_calc_util import get_mds42_rbs_strength_profile
from tool_wrappers.rbs_calc_util import is_rbs_strength_conserved


def check_all(original_seq_record, refactored_seq_record,
        ignore_problems_in_feature_ids=[], interval=None):
    """Verifies obvious rules including translation conservation and RNA
    conservation.

    Does NOT check that all forbidden codons are removed. Use
    check_recoding_successful() for a more extensive set of checks.
    """
    check_translations_conserved(original_seq_record, refactored_seq_record,
            interval=interval)
    check_rnas_conserved(original_seq_record, refactored_seq_record,
            ignore_problems_in_feature_ids=ignore_problems_in_feature_ids,
            interval=interval)

    # NOTE: check_rbs() ommitted because it takes a long time, so clients
    # should call it manually when desired.


def check_recoding_is_complete(original_seq_record, refactored_seq_record,
        ignore_problems_in_feature_ids=[], interval=None):
    """Does more extensive checks, including that all forbidden codons are
    removed.

    Raises:
        AssertionError if anything fails.
    """
    try:
        check_all(original_seq_record, refactored_seq_record,
                ignore_problems_in_feature_ids=ignore_problems_in_feature_ids,
                interval=interval)
        check_forbidden_codons_removed(refactored_seq_record, CODONS_TO_REMOVE)
    except AssertionError as e:
        raise AssertionError(
                "The passed in genome_record has problem: %s" % str(e))


def _get_features_passing_interval_filter(seq_record, interval):
    if interval:
        return [f for f in seq_record.features if
                does_interval_overlap_feature(interval, f)]
    else:
        return seq_record.features


def check_translations_conserved(original_seq_record, refactored_seq_record,
        interval=None):
    """Confirms that the translations of the coding features are preserved
    before and after.
    """
    print '...Checking translation is conserved...'

    original_coding_features = _get_features_passing_interval_filter(
            original_seq_record, interval)
    original_coding_features = filter(
            lambda feature: feature.type == 'CDS',
            original_coding_features)

    refactored_coding_features = _get_features_passing_interval_filter(
            refactored_seq_record, interval)
    refactored_coding_features = filter(
            lambda feature: feature.type == 'CDS',
            refactored_coding_features)

    error_msg = "Different number of CDS features."
    assert len(original_coding_features) == len(refactored_coding_features), error_msg

    for original_feature in original_coding_features:
        if original_feature.id in TRANSLATION_CONSERVED_EXCEPTIONS:
            continue
        original_feature_seq = original_feature.extract(original_seq_record.seq)
        original_translation = translate_custom(str(original_feature_seq))

        refactored_feature = get_feature_by_id(
                refactored_seq_record, original_feature.id)
        if not refactored_feature:
            raise AssertionError("Feature lost after refactor: %s" %
                    str(original_feature))
        refactored_feature_seq = refactored_feature.extract(
                refactored_seq_record.seq)
        try:
            refactored_translation = translate_custom(str(refactored_feature_seq))
        except TranslationError as e:
            print "Error translating %s" % str(original_feature)
            raise e

        error_msg = "Translation mismatch for feature %s" % original_feature.id
        assert original_translation == refactored_translation, error_msg

    print '......Translation conservation confirmed.'

    # All tests passed.
    return True


def check_rnas_conserved(original_seq_record, refactored_seq_record,
        ignore_problems_in_feature_ids=[], interval=None):
    """Check that RNA coding sequences are conserved.
    """
    print '...Checking RNA\'s are conserved...'

    RNA_TYPES = set([
        'misc_RNA',
        'ncRNA',
        'rRNA',
        'tRNA',
        'tmRNA'
    ])

    # Filter original features to those in interval, if provided.
    if interval:
        interval_features = [f for f in original_seq_record.features if
                does_interval_overlap_feature(interval, f)]
    else:
        interval_features = original_seq_record.features
    original_rna_features = filter(
            lambda feature: feature.type in RNA_TYPES and
                    not feature.id in ignore_problems_in_feature_ids,
            interval_features)

    if interval:
        refactored_interval_features = [f for f in
                refactored_seq_record.features if
                does_interval_overlap_feature(interval, f)]
    else:
        refactored_interval_features = refactored_seq_record.features
    refactored_rna_features = filter(
            lambda feature: feature.type in RNA_TYPES and
                    not feature.id in ignore_problems_in_feature_ids,
            refactored_interval_features)

    error_msg = "Different number of RNA features."
    assert len(original_rna_features) == len(refactored_rna_features), error_msg

    for original_feature in original_rna_features:
        original_feature_seq = original_feature.extract(original_seq_record.seq)

        refactored_feature = get_feature_by_id(
                refactored_seq_record, original_feature.id)
        refactored_feature_seq = refactored_feature.extract(
                refactored_seq_record.seq)

        error_msg = "RNA not conserved for %s" % original_feature.id
        assert str(original_feature_seq) == str(refactored_feature_seq), error_msg

    print '......RNA conservation confirmed.'

    return True


def check_rbs_conserved(refactored_seq_record):
    """Confirm that the rbs's for each feature are conserved.
    """
    orig_rbs_expression_profile = get_mds42_rbs_strength_profile()

    refactored_coding_features = filter(
            lambda feature: feature.type == 'CDS',
            refactored_seq_record.features)

    num_features = len(refactored_coding_features)
    for feature_index in range(num_features):
        feature = refactored_coding_features[feature_index]

        print 'Checking RBS conservation for feature %d of num_features %d' % (
                    feature_index + 1, num_features)

        assert is_rbs_strength_conserved(
                orig_rbs_expression_profile,
                feature,
                refactored_seq_record.seq)


def check_forbidden_codons_removed(
        genome_record,
        forbidden_codons,
        feature_start_lower_bound=None,
        feature_start_upper_bound=None):
    """Checks CDS features in the SeqRecord for forbidden codons.

    Accepts optional lower and upper bounds on feature starting locations.
    """
    print '...Checking forbidden codons removed...'

    coding_features = filter(lambda f: f.type == 'CDS', genome_record.features)

    # May apply lower and upper bounds.
    if feature_start_lower_bound:
        coding_features = filter(
                lambda f: f.location.start > feature_start_lower_bound,
                coding_features)
    if feature_start_upper_bound:
        coding_features = filter(
                lambda f: f.location.start < feature_start_upper_bound,
                coding_features)

    # Test for forbidden codons.
    for feature in coding_features:
        # print 'Checking:', feature.id
        num_codons = len(feature) / 3
        seq = str(feature.extract(genome_record.seq)).upper()
        # NOTE: Intentionally skip the starting codon.
        for codon_index in range(1, num_codons):
            codon = seq[codon_index * 3 : codon_index * 3 + 3]

            if SELENOCYSTEINE_CODON == codon:
                # TGA occurring mid-feature, codes for Selenocysteine.
                continue

            assert not codon in forbidden_codons, (
                    "forbidden %s in feature %s at codon index %d" %
                            (codon, feature.id, codon_index)
            )

    print '......No forbidden codons found.'


def check_each_cds_ends_with_stop_codon(seq_record, ref_seq_record,
        stop_codons=['TAA', 'TGA', 'TAG']):
    """Checks whether any stop codons have been unexpectedly lost
    in refactoring a genome.

    Args:
        seq_record: The SeqRecord we are checking.
        ref_seq_record: This original SeqRecord. This is used to exclude cases
            where a particular feature already did not end in a standard stop
            codon.
        stop_codons: Allows overriding the stop codons that we check.

    Returns:
        List of ids of CDS features that don't end with a valid stop codon.
    """
    def _get_feature_last_codon(feature, feature_seq_record):
        feature_seq = feature.extract(feature_seq_record.seq)
        return str(feature_seq[-3:]).upper()

    # Build the reference map.
    ref_seq_gene_to_last_codon_map = {}
    for feature in ref_seq_record.features:
        if feature.type == 'CDS' and 'gene' in feature.qualifiers:
            ref_seq_gene_to_last_codon_map[feature.qualifiers['gene'][0]] = (
                    _get_feature_last_codon(feature, ref_seq_record))

    cds_features_with_invalid_stop_codon = []
    for feature in seq_record.features:
        # Run the gauntlet of checks.
        if not feature.type == 'CDS':
            continue
        last_codon = _get_feature_last_codon(feature, seq_record)
        if last_codon in stop_codons:
            continue
        if ('gene' in feature.qualifiers and
                not ref_seq_gene_to_last_codon_map[feature.qualifiers['gene'][0]]
                        in stop_codons):
                continue

        # Gauntlet passed. Add this feature to the invalid list.
        cds_features_with_invalid_stop_codon.append(feature.id)
    return cds_features_with_invalid_stop_codon


def check_feature_for_internal_stops(feature, seq_record,
        stop_codons=set(['TAA', 'TGA', 'TAG'])):
    """Checks a single feature for internal stop codons.

    Args:
        feature: A SeqFeature.
        seq_record: SeqRecord that feature comes from. Used to get the
            feature sequence.

    Returns:
        A list of codon indeces where the feature has internal stops.
    """
    feature_bad_codons = []
    num_codons = len(feature) / 3
    seq = str(feature.extract(seq_record.seq))
    # Check all except last codon.
    for codon_index in range(1, num_codons - 1):
        codon = seq[codon_index * 3 : codon_index * 3 + 3]
        if codon in stop_codons:
            feature_bad_codons.append(str((codon, codon_index)))
    return feature_bad_codons


def check_internal_stops(seq_record, report_file,
        stop_codons=set(['TAA', 'TGA', 'TAG'])):
    """Checks the features of seq_record for internal stop codons.
    """
    num_with_internal_stops = 0
    with open(report_file, 'w') as output_fh:
        for feature in seq_record.features:
            # Run the gauntlet of checks.
            if not feature.type == 'CDS':
                continue

            if 'pseudo' in feature.qualifiers:
                continue

            feature_bad_codons = check_feature_for_internal_stops(
                    feature, seq_record, stop_codons=stop_codons)

            # Print results
            if len(feature_bad_codons) > 0:
                num_with_internal_stops += 1
                if 'gene' in feature.qualifiers:
                    print feature.qualifiers['gene'][0]
                else:
                    print feature
                output_fh.write('Feature: %s\n, internal stops: %d, %s' % (
                        str(feature), len(feature_bad_codons),
                        str(feature_bad_codons)))

    print 'Total: ', num_with_internal_stops


def check_cds_start_codons(seq_record):
    """Checks the features of seq_record for internal stop codons.
    """
    for feature in seq_record.features:
        # Run the gauntlet of checks.
        if not feature.type == 'CDS':
            continue

        if 'pseudo' in feature.qualifiers:
            continue

        if not does_feature_start_with_valid_start_codon(
                feature, seq_record):
            print feature


if __name__ == '__main__':
    import os

    import yaml

    from biopython_util import get_genome_record
    from paths import CONFIG_DIR
    from paths import OUTPUT_DIR

    from refactor_config import ORIGINAL_GENOME_RECORD

    IGNORE_PROBLEMS_IN_FEATURE_IDS = [
            # Homopolymers
            'ECMDS42_0735',
            'ECMDS42_0867',

            # Restriction sites
            'ECMDS42_1514',
            'ECMDS42_1602',
            'ECMDS42_2065'
    ]

    MDS42_RECORD = ORIGINAL_GENOME_RECORD


    CONFIG = os.path.join(CONFIG_DIR, 'seven_codon_config.yaml')
    with open(CONFIG) as yaml_fh:
        YAML_CONFIG_DICT = yaml.load(yaml_fh)
    CODONS_TO_REMOVE = YAML_CONFIG_DICT['forbidden_codons']

    REFACTORED_GENOME = '/home/glebk/Projects/churchlab/genome-refactor/data/completed_segments/seg5-seg48/2013_08_30.seg5-seg48.motifs_removed.AGN_fixed.muddled.genbank'

    REFACTORED_GENOME_RECORD = get_genome_record(REFACTORED_GENOME)

    check_translations_conserved(MDS42_RECORD, REFACTORED_GENOME_RECORD)
    check_rnas_conserved(MDS42_RECORD, REFACTORED_GENOME_RECORD, ignore_problems_in_feature_ids=IGNORE_PROBLEMS_IN_FEATURE_IDS)

    check_forbidden_codons_removed(
            REFACTORED_GENOME_RECORD,
            CODONS_TO_REMOVE)
