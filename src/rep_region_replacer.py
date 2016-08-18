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
Methods for replacing REP regions with terminators.
"""

from collections import defaultdict
import copy
import csv
import os
import re

from Bio.Seq import reverse_complement
from Bio.SeqFeature import SeqFeature
import pandas as pd

from biopython_util import add_feature_to_seq_record
from biopython_util import calc_interval_list_to_features_overlapped
from biopython_util import delete_interval
from biopython_util import does_interval_overlap_feature
from biopython_util import insert_sequence_and_update_features
from biopython_util import get_feature_label
from biopython_util import set_feature_label
from paths import DATA_DIR
from post_processing.gen9_synthesis_constraints_checker import check_all



# Annotations
REP_REGION_ANNOTATION = 'repeat_region'
MARKED_REP_REGION_FOR_DELETE = 'rep_to_delete'
REPLACEMENT_TERMINATOR = 'rep_to_term'

# Source data
VOIGT_TERMINATOR_DATA = os.path.join(DATA_DIR,
        'voigt_2013_synthetic_terminators.csv')

# Filtered - see generate_filtered_terminators()
VOIGT_TERMINATOR_DATA_FILTERED = os.path.join(DATA_DIR,
        'voigt_2013_synthetic_terminators_filtered.csv')


def mark_rep_regions_to_delete(genome_record, interval_list, ignore_rep_list,
        min_rep_size=50):
    """Returns modified genome_record with rep regions to delete marked.

    Args:
        genome_record: SeqRecord with features.
        interval_list: Only repeats overlapping these intervals will be
            marked.
        ignore_rep_list: List of REP numbers to ignore.
    """
    ignore_rep_set = set(ignore_rep_list)

    # Only keep REPs
    repeat_region_features = [feature for feature in genome_record.features
            if feature.type == REP_REGION_ANNOTATION]

    # Only keep those in valid regions.
    relevant_repeat_regions = []
    for repeat_region_feature in repeat_region_features:
        if len(repeat_region_feature) < min_rep_size:
            continue
        for interval in interval_list:
            is_overlap = does_interval_overlap_feature(interval,
                    repeat_region_feature)
            if (is_overlap and
                    not _get_rep_number(repeat_region_feature)
                            in ignore_rep_set):
                relevant_repeat_regions.append(repeat_region_feature)

    # Grab cds features for checking overlaps below.
    cds_features = [f for f in genome_record.features if f.type == 'CDS']
    assert len(cds_features) > 0

    # First we mark the features annotated repeat regions.
    for repeat_region_feature in relevant_repeat_regions:
        feature_id = _get_rep_id(repeat_region_feature)

        # Check that the rep region doesn't overlap any cds features.
        rep_interval = (
            repeat_region_feature.location.start,
            repeat_region_feature.location.end)
        all_overlaps = calc_interval_list_to_features_overlapped(
                [rep_interval], cds_features)
        for o in all_overlaps:
            if len(o):
                raise AssertionError("REP %d overlaps %s" % (
                    _get_rep_number(repeat_region_feature), str(o)))

        rep_region_to_delete = SeqFeature(
                location=copy.deepcopy(repeat_region_feature.location),
                type=MARKED_REP_REGION_FOR_DELETE,
                strand=1,
                id=feature_id
        )
        set_feature_label(rep_region_to_delete, feature_id)
        add_feature_to_seq_record(genome_record, rep_region_to_delete)


def _get_rep_id(feature):
    feature_note = feature.qualifiers['note'][0]
    return 'd_' + re.match(r'(R[E|I]P[0-9]+)', feature_note).group()


def _get_rep_number(feature):
    feature_note = feature.qualifiers['note'][0]
    return int(re.match(r'R[E|I]P([0-9]+)', feature_note).group(1))


def UniqueSyntheticTerminatorGenerator(max_terminator_strength_ts=200):
    """Yields a DataFrame series with data about a particular terminator.

    Yields from strongest to weakest as determined by order of the input csv.

    Raises:
        StopIteration when all terminators have been checked.
    """
    voigt_synthetic_terminator_df = pd.read_csv(VOIGT_TERMINATOR_DATA)
    index = 0
    while index < len(voigt_synthetic_terminator_df):
        candidate = voigt_synthetic_terminator_df.iloc[index]
        index += 1
        if candidate['Average Strength'] > max_terminator_strength_ts:
            continue

        # Verify against Gen9 constraints.
        seq_verification_result = check_all(candidate['Sequence'])
        verified = True
        for key in seq_verification_result.iterkeys():
            if len(seq_verification_result[key]):
                verified = False
                break
        if not verified:
            continue

        yield candidate


def FilteredCyclicSyntheticTerminatorGenerator(first_terminator_seq=None):
    """Yields a DataFrame series with data for a terminator from the filtered
    set, cycling from the beginning after all terminators have been used.

    Args:
        first_terminator_seq: If provided, skip first n terminators that are
            not this one, and then resume iterating as normal.
    """
    filtered_terminators_df = pd.read_csv(VOIGT_TERMINATOR_DATA_FILTERED)
    num_terminators = len(filtered_terminators_df)
    index = 0
    yielded_at_least_one = False
    while True:
        # Check whether we need to skip through to the first iterator.
        if first_terminator_seq is not None and not yielded_at_least_one:
            while not yielded_at_least_one:
                assert index < num_terminators, (
                        "Invalid first_terminator_seq" % first_terminator_seq)
                if (filtered_terminators_df.iloc[index]['Sequence'] ==
                        first_terminator_seq):
                    break
                index += 1
        yield filtered_terminators_df.iloc[index % num_terminators]
        yielded_at_least_one = True
        index += 1


def assign_terminators(rep_to_terminator_df,
        terminator_generator=UniqueSyntheticTerminatorGenerator()):

    def _assign_terminators(row):
        if row['terminator_dir'] == 0:
            return pd.Series({})

        t_data = terminator_generator.next()
        if abs(row['terminator_dir']) > 0:
            result_dict = {
                'upstream_terminator_name': t_data['Name'],
                'upstream_terminator_strength_avg': t_data[
                        'Average Strength'],
                'upstream_terminator_sequence': t_data[
                        'Sequence']
            }
            if row['terminator_dir'] == 2:
                t_data = terminator_generator.next()
                result_dict.update({
                    'downstream_terminator_name': t_data['Name'],
                    'downstream_terminator_strength_avg': t_data[
                            'Average Strength'],
                    'downstream_terminator_sequence': t_data[
                            'Sequence']
                })

        return pd.Series(result_dict)

    return rep_to_terminator_df.join(
            rep_to_terminator_df.apply(_assign_terminators, axis=1))


def delete_marked_rep_regions(genome_record, rep_to_terminator_df):
    # TODO: What's the pandas way to do the following?
    # I think I did it up above with set_index(). Had some int vs string
    # snafoo here.
    rep_id_to_terminator_map = {}
    for idx, row in rep_to_terminator_df.iterrows():
        rep_id_to_terminator_map[row['rep_id']] = row

    next_rep_to_delete = _get_next_rep_marked_for_deletion(genome_record)
    while (next_rep_to_delete):
        delete_rep_id = int(get_rep_id_from_delete_rep_id(
                get_feature_label(next_rep_to_delete)))
        print 'deleting %s' % delete_rep_id

        rep_interval = (next_rep_to_delete.location.start,
                next_rep_to_delete.location.end)
        delete_interval(genome_record, rep_interval)

        # Maybe replace with terminator.
        replacement_terminator_data = rep_id_to_terminator_map[delete_rep_id]
        if abs(replacement_terminator_data['terminator_dir']) >= 1:
            # Add upstream (or only) terminator.
            upstream_terminator_seq = replacement_terminator_data[
                    'upstream_terminator_sequence']
            strand = 1
            if replacement_terminator_data['terminator_dir'] == -1:
                upstream_terminator_seq = reverse_complement(
                        upstream_terminator_seq)
                strand = -1
            feature_id = str(delete_rep_id) + '_upstream_terminator'
            insert_sequence_and_update_features(
                    genome_record,
                    upstream_terminator_seq,
                    rep_interval[0],
                    insert_feature_type=REPLACEMENT_TERMINATOR,
                    insert_feature_id=feature_id,
                    insert_feature_strand=strand
            )
            if replacement_terminator_data['terminator_dir'] == 2:
                # Add downstream reverse terminator.
                downstream_terminator_seq = replacement_terminator_data[
                        'downstream_terminator_sequence']
                downstream_terminator_seq = reverse_complement(
                        downstream_terminator_seq)
                strand = -1
                feature_id = str(delete_rep_id) + '_downstream_terminator'
                insert_sequence_and_update_features(
                        genome_record,
                        downstream_terminator_seq,
                        rep_interval[0] + len(upstream_terminator_seq),
                        insert_feature_type=REPLACEMENT_TERMINATOR,
                        insert_feature_id=feature_id,
                        insert_feature_strand=strand
                )

        # Next iteration.
        next_rep_to_delete = _get_next_rep_marked_for_deletion(genome_record)


def _get_next_rep_marked_for_deletion(genome_record):
    for feature in genome_record.features:
        if feature.type == MARKED_REP_REGION_FOR_DELETE:
            return feature
    return None


def get_rep_id_from_delete_rep_id(delete_rep_id):
    maybe_match = re.match(r'd_R[E|I]P([0-9]+)', delete_rep_id)
    if maybe_match:
        return maybe_match.group(1)
    return None


def debug_generate_valid_terminators():
    """Generate a list of valid terminators according to above constraints.
    """
    terminator_generator = UniqueSyntheticTerminatorGenerator(
            max_terminator_strength_ts=999999)

    VOIGT_FILTERED_VARIANTS = (
        'voigt_terminators_filtered_by_synthesis_constraints.csv')
    FIELD_NAMES = [
        'name',
        'strength_avg',
        'strength_std',
        'sequence',
        'A_tract',
        'hairpin',
        'structure',
        'U_tract',
    ]

    with open(VOIGT_FILTERED_VARIANTS, 'w') as csv_fh:
        writer = csv.DictWriter(csv_fh, FIELD_NAMES)
        writer.writeheader()
        for t_data in terminator_generator:
            writer.writerow({
                'name': t_data['Name'],
                'strength_avg': t_data['Average Strength'],
                'sequence': t_data['Sequence'],
                'A_tract': t_data['A-tract'],
                'hairpin': t_data['Hairpin'],
                'structure': t_data['Structure'],
                'U_tract': t_data['U-tract'],
            })


def generate_filtered_terminators(min_terminator_strength_ts=5):
    """Offline function that generates a filtered set of terminators.

    Filters include:
        * Must pass synthesis constraints
        * Minimum strength
        * Get rid of repeat hairpin sequences (choose one with average) strength

    Terminator sequences adjusted to only be the A-tract, hairpin, and U-tract.
    """
    voigt_synthetic_terminator_df = pd.read_csv(VOIGT_TERMINATOR_DATA)

    # First, we modify the terminator sequences from what the Voigt data has.
    # Specifically, the 'Sequence' column in the Voigt data has a few issues:
    #     * Extra flanking context sequences on some constructs that is reused
    #       in many of the constructs, which would lead to additional homology
    #       if we are to use it for replacing REP regions.
    #     * Sometimes the Sequence column is missing the A-tract or U-tract
    # We make the assumption that we can remove the context without affecting
    # the strength too much, and then combine the A-tract, hairpin, and U-tract
    # which gets rid of the buggy 'Sequence' column issue.
    voigt_synthetic_terminator_df['Sequence'] = (
            voigt_synthetic_terminator_df['A-tract'] +
            voigt_synthetic_terminator_df['Hairpin'] +
            voigt_synthetic_terminator_df['U-tract'])
    voigt_synthetic_terminator_df['Sequence'] = (
            voigt_synthetic_terminator_df['Sequence'].apply(
                    lambda x: x.replace('U', 'T')))

    # Filter out those where average strength is too low.
    voigt_synthetic_terminator_df = voigt_synthetic_terminator_df[
            voigt_synthetic_terminator_df['Average Strength'] >
                    min_terminator_strength_ts]

    # Filter out terminator sequences that don't pass synthesis constraints.
    def _check_synthesis_constraints(sequence):
        """Checks whether the terminator sequence passes synthesis constraints.
        """
        seq_verification_result = check_all(sequence)
        verified = True
        for key in seq_verification_result.iterkeys():
            if len(seq_verification_result[key]):
                verified = False
                break
        return verified
    voigt_synthetic_terminator_df = voigt_synthetic_terminator_df[
        voigt_synthetic_terminator_df['Sequence'].apply(
                _check_synthesis_constraints)]

    # Since hairpins are repeated, select only one row for each unique hairpin,
    # using the one with the median strength.
    hairpin_to_terminator_data_map = defaultdict(list)
    for idx, row in voigt_synthetic_terminator_df.iterrows():
        hairpin_to_terminator_data_map[row['Hairpin']].append(row)
    unique_hairpin_rows = []
    for hairpin, rows in hairpin_to_terminator_data_map.iteritems():
        # Sort in decreasing order so we always select the one having the
        # median, or of the next greater strength.
        sorted_rows = sorted(rows, key=lambda x: x['Average Strength'],
                reverse=True)
        middle_idx = len(sorted_rows) / 2
        unique_hairpin_rows.append(sorted_rows[middle_idx])
    voigt_synthetic_terminator_df = pd.DataFrame(unique_hairpin_rows)

    # Write results to output file.
    voigt_synthetic_terminator_df = voigt_synthetic_terminator_df[[
        'Name',
        'Average Strength',
        'Sequence',
        'A-tract',
        'Hairpin',
        'Structure',
        'U-tract'
    ]]
    voigt_synthetic_terminator_df.to_csv(VOIGT_TERMINATOR_DATA_FILTERED,
            index=False)


# Map to deal with our input file.
FLANKING_GENE_ORIENTATION_TO_TERMINATOR_DIR = {
    ('forward', 'forward'): 1,
    ('reverse', 'reverse'): -1,
    ('forward', 'reverse'): 2, # terminators in both directions. t>>>> <<<<t
    ('reverse', 'forward'): 0,
}


def generate_rep_id_to_replacement_terminator_orientation_map(
        genome_record, rep_gene_orientation_data):
    """Returns Dataframe mapping rep id to the orientation that the
    replacement terminator should have.
    """
    # First find the list of REPs we are deleting.
    marked_for_delete_features = [feature for feature in genome_record.features
            if feature.type == MARKED_REP_REGION_FOR_DELETE]

    rep_data_list = [{
        'rep_id': get_rep_id_from_delete_rep_id(get_feature_label(feature)),
        'rep_size': len(feature),
    } for feature in marked_for_delete_features]
    deleted_rep_id_df = pd.DataFrame(rep_data_list)

    # Data about orientation of genes flanking rep regions.
    rep_orientation_df = pd.read_csv(rep_gene_orientation_data)
    rep_orientation_df = rep_orientation_df.set_index('rep_id')

    # We build the mapping using the following helper function.
    def _add_terminator_data(row):
        """Helper function that determines orientation for replacement terminator
        depending on flanking genes.
        """
        # Figure out id.
        rep_id = row['rep_id']
        rep_size = row['rep_size']

        # Determine replacement terminator orientation.
        upstream_dir = rep_orientation_df.loc[rep_id][
                'upstream_gene_orientation']
        downstream_dir = rep_orientation_df.loc[rep_id][
                'downstream_gene_orientation']
        terminator_dir = FLANKING_GENE_ORIENTATION_TO_TERMINATOR_DIR[
                (upstream_dir, downstream_dir)]

        if terminator_dir == 2:
            # Only do the direction of the essential gene, unless they
            # are both essential.
            upstream_essentiality = rep_orientation_df.loc[rep_id][
                'upstream_gene_essentiality']
            downstream_essentiality = rep_orientation_df.loc[rep_id][
                'downstream_gene_essentiality']
            if upstream_essentiality == 'E' and downstream_essentiality == 'E':
                terminator_dir = 2
            elif upstream_essentiality == 'E':
                terminator_dir = 1
            elif downstream_essentiality == 'E':
                terminator_dir = -1
            else:
                # Forward-facing terminator by default.
                # TODO: Or do we want no terminator.
                terminator_dir = 1

        return {
            'rep_id': rep_id,
            'rep_size': rep_size,
            'terminator_dir': terminator_dir
        }

    # NOTE: Rewrote as a iterrows because I couldn't debug when this was
    # executed as an apply() call.
    mapping_data_list = []
    for idx, row in deleted_rep_id_df.iterrows():
        mapping_data_list.append(_add_terminator_data(row))
    return pd.DataFrame(mapping_data_list)


if __name__ == '__main__':
    generate_filtered_terminators()
