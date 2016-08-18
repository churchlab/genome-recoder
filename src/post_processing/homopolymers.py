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
Methods for dealing with occurrences of homopolymers.

The main constraint in this space we have right now is imposed by Gen9
synthesis requirements where there should not be homopolymer runs of
greater than 5 bases.
"""

import copy
import csv
import re

from biopython_util import COMMONLY_IGNORED_FEATURE_TYPES
from biopython_util import calc_interval_list_to_features_overlapped
from biopython_util import get_region_codon_indeces_in_feature
from biopython_util import update_seq_record_feature
from codon_replacer_util import replace_codons_in_single_feature
from feature_profile import SecondaryStructureFeatureProfile
from refactor_config import CODONS_TO_REMOVE
from refactor_config import HOMOPOLYMER_RUN_LIMIT_A
from refactor_config import HOMOPOLYMER_RUN_LIMIT_C
from refactor_config import HOMOPOLYMER_RUN_LIMIT_G
from refactor_config import HOMOPOLYMER_RUN_LIMIT_T
from refactor_config import ORIGINAL_CODON_USAGE_MEMEX
from refactor_config import REFACTORED_CODON_USAGE_MEMEX


# Regular expression used for finding runs of homopolymers.
HOMOPOLYMER_RE = (
        '(?P<A>A{' + str(HOMOPOLYMER_RUN_LIMIT_A + 1) + ',})|'
        '(?P<C>C{' + str(HOMOPOLYMER_RUN_LIMIT_C + 1) + ',})|'
        '(?P<G>G{' + str(HOMOPOLYMER_RUN_LIMIT_G + 1) + ',})|'
        '(?P<T>T{' + str(HOMOPOLYMER_RUN_LIMIT_T + 1) + ',})')


def find_homopolymer_runs(genome_record, start_bound=None, end_bound=None):
    """Finds runs of homopolymer runs, sorted by start.

    Args:
        genome_record: The SeqRecord to search over.
        start_bound: Optionally bound the search to start at this position.
        end_bound: Optionally bound the search to end at this position.

    Returns:
        List of objects with keys:
            * interval: Tuple-pair of (start, end).
            * group: The actual matched run (e.g. 'AAAAAAA').
    """
    seq = str(genome_record.seq)

    # Parse the bounds.
    effective_start_bound = start_bound if start_bound else 0
    effective_end_bound = end_bound if end_bound else len(seq)
    effective_seq = seq[effective_start_bound:effective_end_bound]

    # Search for the results.
    h_run_obj_list = [
            {
                'interval': (i.start() + effective_start_bound,
                        i.end() + effective_start_bound),
                'group': i.group()
            }
            for i in re.finditer(HOMOPOLYMER_RE, effective_seq)]

    # Make sure the results are sorted by start position.
    h_run_obj_list = sorted(h_run_obj_list, key=lambda x: x['interval'][0])

    return h_run_obj_list


def remove_homopolymer_runs(refactor_context, start_bound=None, end_bound=None,
        report_prefix=None):
    """Scans the genome_record and removes instances of homopolyme runs.

    The strategy varies according to whether the homopolymer run is inside of
    a coding feature or not.

    Args:
        refactor_context: The RefactorContext.
        start_bound: Optionally bound fixes to start at this position.
        end_bound: Optionally bound fixes to end at this position.
        report_fix: Destination to write flagged sites.

    Returns:
        Object with the following keys:
            * updated_genome_record: A copy of the genome_record in the
            passed-in refactor_context
            * flagged: List of flagged homopolymer runs.
    """
    print '...Removing homopolymer runs...'

    # Mutable values to be updated and returned.
    mutable_genome_record = copy.deepcopy(refactor_context.get_genome_record())
    flagged = []

    # Identify homopolymer run occurrences.
    homopolymer_runs = find_homopolymer_runs(mutable_genome_record,
            start_bound=start_bound, end_bound=end_bound)
    interval_list = [h['interval'] for h in homopolymer_runs]

    # Identify all overlaps in advance. This allows us to find the overlaps
    # more efficiently, as well as check for multiple overlaps.
    relevant_features = [
            feature for feature in mutable_genome_record.features
            if not feature.type in COMMONLY_IGNORED_FEATURE_TYPES]
    h_run_idx_to_feature_overlapped_list = (
            calc_interval_list_to_features_overlapped(
                    interval_list, relevant_features))

    # Loop over the homopolyer run occurrences and either fix or flag them.
    num_runs = len(homopolymer_runs)
    for h_run_idx, h_run_obj in enumerate(homopolymer_runs):
        print 'Removing homopolymer %d of %d' % (h_run_idx + 1, num_runs)

        features_overlapping_h_run = h_run_idx_to_feature_overlapped_list[
                h_run_idx]

        if len(features_overlapping_h_run) == 0:
            h_run_obj['exception_string'] = 'not annotated'
            flagged.append(h_run_obj)
        elif len(features_overlapping_h_run) >= 2:
            h_run_obj['exception_string'] = 'overlaps more than 1 feature'
            flagged.append(h_run_obj)
        else:
            feature = features_overlapping_h_run[0]
            if feature.type == 'CDS':
                remove_result = (
                        _remove_homopolymer_run_in_coding_feature(
                                refactor_context,
                                mutable_genome_record,
                                h_run_obj,
                                feature))
                if remove_result['is_success']:
                    mutable_genome_record = remove_result[
                            'updated_genome_record']
                else:
                    h_run_obj['exception_string'] = remove_result[
                            'exception_string']
                    flagged.append(h_run_obj)
            else:
                h_run_obj['exception_string'] = (
                        'in RNA-like feature %s' % feature.id)
                h_run_obj['pos_in_feature'] = (h_run_obj['interval'][0] -
                        feature.location.start + 1)
                flagged.append(h_run_obj)

    if report_prefix:
        report_file_output = report_prefix + 'homopolymers.csv'
        REPORT_FIELDNAMES = [
            'exception_string',
            'interval',
            'group',
            'pos_in_feature'
        ]

        with open(report_file_output, 'w') as report_fh:
            writer = csv.DictWriter(report_fh, REPORT_FIELDNAMES)
            writer.writeheader()

            for flagged_site in flagged:
                writer.writerow(flagged_site)


    print '...Done removing homopolymer runs.'
    return {
            'updated_genome_record': mutable_genome_record,
            'flagged': flagged
    }


def _remove_homopolymer_run_in_coding_feature(
        refactor_context,
        mutable_genome_record,
        h_run_obj,
        feature):
    """Removes a homopolymer run in a coding feature.

    The strategy is to muddle up all affected codons within the feature
    so as to reduce the chance of "snap-back" over generations.

    Args:
        refactor_context: The RefactorContext.
        mutable_genome_record: The SeqRecord object representing the genome.
        h_run_obj: Object containing data about a specific occurrence of a
            homopolyer run.
        feature: The feature that is overlapped by the restriction site.

    Returns:
        Object with keys:
            * is_success: Whether remove succeeded.
            * updated_genome_record: The updated genome_record if successful.
            * exception_string: Message describing failure.
    """
    # First check if this is only a partial overlap and thus we can avoid
    # a heavy change.
    interval = h_run_obj['interval']
    interval_start = interval[0]
    interval_end = interval[1]
    interval_size = interval_end - interval_start

    # Otherwise commence the muddling strategy.
    feature_seq = str(feature.extract(mutable_genome_record.seq))

    # Figure out the specific codons that need to be changed.
    affected_codon_indeces = get_region_codon_indeces_in_feature(
            feature, interval)
    avoid_codons_in_positions = {}
    for codon_index in affected_codon_indeces:
        codon = feature_seq[codon_index * 3 : codon_index * 3 + 3]
        avoid_codons_in_positions[codon_index] = codon

    # Perform replace.
    first_codon_to_modify = affected_codon_indeces[0]
    last_codon_to_modify = affected_codon_indeces[-1]
    assert first_codon_to_modify <= last_codon_to_modify
    result = replace_codons_in_single_feature(
            refactor_context,
            feature.id,
            explicit_genome_record=mutable_genome_record,
            start_codon_index=first_codon_to_modify,
            last_codon_index=last_codon_to_modify,
            avoid_codons_in_positions=avoid_codons_in_positions)
    if not result['is_success']:
        return {
                'is_success': False,
                'exception_string': result['exception_string']
        }

    update_seq_record_feature(
            mutable_genome_record,
            feature.id,
            result
    )
    return {
            'is_success': True,
            'updated_genome_record': mutable_genome_record
    }


if __name__ == '__main__':
    import cProfile
    from datetime import datetime
    import os
    import pickle

    from Bio import SeqIO

    from biopython_util import get_genome_record
    from refactor_config import OUTPUT_DIR
    from refactor_context import RefactorContext

    genome_record = get_genome_record(
            '../data/2013_03_06_20_16_04_mds42_refactored.gbk')
    START_BOUND = 100000
    END_BOUND = 150000
    # runs = find_homopolymer_runs(genome_record, start_bound=START_BOUND,
    #         end_bound=END_BOUND)
    # print runs

    refactor_context = RefactorContext(genome_record)

    print 'Removing homopolymers...'
    tmp_file_prefix = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    cprofile_output_dest = os.path.join(
            OUTPUT_DIR, tmp_file_prefix + 'remove_homopolymers_cprofile.out')
    cmd = ('remove_homopolymer_result = remove_homopolymer_runs('
            'refactor_context, start_bound=START_BOUND, end_bound=END_BOUND)'
    )
    cProfile.run(cmd, cprofile_output_dest)

    updated_genome_record = remove_homopolymer_result['updated_genome_record']
    flagged_h_runs = remove_homopolymer_result['flagged']

    # Write the resulting genome.
    genome_output_file = (
            '../data/2013_03_06_20_16_04_mds42_refactored_homopolymers_removed.gbk'
    )
    with open(genome_output_file, 'w') as output_fh:
        SeqIO.write(updated_genome_record, output_fh, 'genbank')

    # Write the flagged data.
    metadata_output_file = (
            '../data/2013_03_06_20_16_04_mds42_refactored_homopolymers_removed.flagged'
    )
    with open(metadata_output_file, 'w') as metadata_output_fh:
        pickle.dump(flagged_h_runs, metadata_output_fh)
