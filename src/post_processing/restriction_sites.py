"""
Methods for finding/removing restriction sites.
"""

import copy
import csv

from Bio.Seq import Seq
from Bio.Restriction import FormattedSeq

from biopython_util import COMMONLY_IGNORED_FEATURE_TYPES
from biopython_util import does_interval_overlap_feature
from biopython_util import get_region_codon_indeces_in_feature
from biopython_util import update_seq_record_feature
from codon_replacer_util import replace_codons_in_single_feature
from refactor_config import RESTRICTION_ENZYME_SITES_TO_REMOVE


def find_restriction_site_occurrences(genome_record, enzyme,
        start_bound=None, end_bound=None):
    """Finds each occurrence of a restriction site across the genome_record.

    Args:
        genome_record: A SeqRecord object representing the genome.
        enzyme: A Bio.Restriction.RestrictionType object representing a
            particular site (e.g. BsmBI).

    Returns:
        A list of objects containing information about where a specific site
        occurs and which strand it occurs on, with keys:
            * interval: Pythonic (start, end) of the site relative to the
                forward strand of genome_record.seq.  The first base occurs
                at start and the last base at end - 1.
            * strand: Whether the the site occurs in the forward or reverse
                reading frame.

    """
    # Parse the bounds.
    seq = genome_record.seq
    effective_start_bound = start_bound if start_bound else 0
    effective_end_bound = end_bound if end_bound else len(seq)
    effective_seq = seq[effective_start_bound:effective_end_bound]

    occurrences = find_restriction_sites_in_seq(effective_seq, enzyme)

    # Update occurences with effective bounds.
    for occur in occurrences:
        orig_interval = occur['interval']
        occur['interval'] = (orig_interval[0] + effective_start_bound,
                orig_interval[1] + effective_start_bound)

    return occurrences


def find_restriction_sites_in_seq(seq, enzyme):
    occurrences = []

    # Search strategy similar to Bio.Restriction.RestrictionBatch.search(),
    # which is to wrap the sequence we're searching across FormattedSeq object
    # and then use the finditer() method.
    forward_group_name = str(enzyme)
    if not isinstance(seq, Seq):
        seq = Seq(seq)
    fseq = FormattedSeq(seq)
    occurrence_iterator = fseq.finditer(enzyme.compsite, enzyme.size)

    occurrence_iterator = fseq.finditer(enzyme.compsite, enzyme.size)
    for start, match_group in occurrence_iterator:
        interval = (start, start + enzyme.size)
        if match_group(forward_group_name):
            strand = 1
        else:
            strand = -1

        # Adjust the interval to be pythonic.
        interval = tuple([pos - 1 for pos in interval])

        # Append the data object representing this occurrence.
        site_occur_data = {
                'enzyme': forward_group_name,
                'site': enzyme.site,
                'interval': interval,
                'strand': strand
        }
        occurrences.append(site_occur_data)

    return occurrences


def find_all_restriction_sites_in_seq(seq, restriction_sites_of_interest):
    all_res_site_occurrences = []
    for enzyme in restriction_sites_of_interest:
        all_res_site_occurrences.extend(find_restriction_sites_in_seq(
                seq, enzyme))
    return all_res_site_occurrences


def remove_restriction_sites(refactor_context, restriction_sites_to_remove,
        start_bound=None, end_bound=None, report_prefix=None):
    """Scans through the entire sequence and removes restriction sites.

    If the site falls inside of an annotated feature, perform synonymous
    codon swapping. Otherwise just muddle it up.

    Nuances:
        * Site may lie in the forward or reverse frame.

    Strategy:
        for each restiction site:
            - find all (start, end) intervals where the site occurs.
            - if the interval falls in an annotated region:
                - if the region is a CDS:
                    - perform synonymous swapping over the overlapping portion.
                - elif the region is an RNA:
                    - if any of the site lies outside of the region:
                        - toggle using rule: A <-> T, or G <-> C

    Args:
        refactor_context: The RefactorContext.
        restriction_sites_to_remove: List of Bio.Restriction.RestrictionType
            objects who site should be removed from the GenomeRecord.
        start_bound: Optional upstream bound for sites to be removed.
        end_bound: Optional downstream bound for sites to be removed.

    Returns:
        Object with keys:
            * updated_genome_record: A copy of the passed in genome_record,
                but with restriction sites removed.
            * flagged: Dictionary from site name to list of flagged
                instances that need to be manually dealt with.
    """
    print '...Removing restriction enzymes...'
    mutable_genome_record = copy.deepcopy(refactor_context.get_genome_record())
    flagged = {}
    for enzyme in restriction_sites_to_remove:
        result = remove_enzyme_site(
                refactor_context, mutable_genome_record, enzyme,
                start_bound=start_bound, end_bound=end_bound,
                report_prefix=report_prefix)
        mutable_genome_record = result['updated_genome_record']
        flagged[str(enzyme)] = result['flagged']
    print '...Done removing restriction enzymes.'
    return {
            'updated_genome_record': mutable_genome_record,
            'flagged': flagged
    }


def remove_enzyme_site(refactor_context, mutable_genome_record, enzyme,
        start_bound=None, end_bound=None, report_prefix=None):
    """Removes every instance of the site associated with the enzyme
    from the genome.

    If the site falls inside of an annotated feature, perform synonymous
    codon swapping. Otherwise just muddle it up. If it's partially inside
    of a feature, change a base that's outside of it.

    Returns:
        The arg mutable_genome_record which was likely mutated.
    """
    flagged = []

    relevant_features = [
            feature for feature in mutable_genome_record.features
            if not feature.type in COMMONLY_IGNORED_FEATURE_TYPES]

    site_occurrences = find_restriction_site_occurrences(
            mutable_genome_record, enzyme, start_bound=start_bound,
            end_bound=end_bound)
    num_occurrences = len(site_occurrences)
    for site_occur_obj_idx in range(num_occurrences):
        print 'Removing %s occurrence %d of %d' % (
                str(enzyme), site_occur_obj_idx + 1, num_occurrences)
        site_occur_obj = site_occurrences[site_occur_obj_idx]
        interval = site_occur_obj['interval']
        at_least_one_overlap_found = False
        for feature in relevant_features:
            is_overlap = does_interval_overlap_feature(interval, feature)
            if is_overlap:
                if feature.type == 'CDS':
                    remove_result = _remove_site_in_coding_feature(
                            refactor_context,
                            mutable_genome_record,
                            site_occur_obj,
                            feature)
                    if remove_result['is_success']:
                        mutable_genome_record = remove_result[
                                'updated_genome_record']
                    else:
                        site_occur_obj['exception_string'] = remove_result[
                                'exception_string']
                        flagged.append(site_occur_obj)
                else:
                    # site_occur_obj['exception_string'] = (
                    #         'in RNA-like feature %s, %s' % (
                    #                 feature.id, feature.qualifiers['gene'][0]))
                    site_occur_obj['exception_string'] = (
                            'in RNA-like feature %s' % feature.id)
                    site_occur_obj['pos_in_feature'] = (
                            site_occur_obj['interval'][0] -
                            feature.location.start + 1)
                    flagged.append(site_occur_obj)
                at_least_one_overlap_found = True
                break

        # Otherwise it is not in an annotated feature.
        if not at_least_one_overlap_found:
            site_occur_obj['exception_string'] = 'not annotated'
            flagged.append(site_occur_obj)

    if report_prefix:
        report_file_output = report_prefix + 'res_site.' + str(enzyme) + '.csv'
        REPORT_FIELDNAMES = [
            'exception_string',
            'interval',
            'strand',
            'pos_in_feature',
            'site',
            'enzyme'
        ]

        with open(report_file_output, 'w') as report_fh:
            writer = csv.DictWriter(report_fh, REPORT_FIELDNAMES)
            writer.writeheader()

            for flagged_site in flagged:
                writer.writerow(flagged_site)

    return {
            'updated_genome_record': mutable_genome_record,
            'flagged': flagged
    }


NON_CODING_TRANSITION_TABLE = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
}


def _remove_site_not_in_coding_feature(mutable_genome_record, site_occur_obj):
    """Removes a restriction site that falls in a region that is not
    annotated as coding.

    The strategy is just to toggle the first base.
    """
    interval = site_occur_obj['interval']
    interval_start = interval[0]
    orig_seq = mutable_genome_record.seq
    updated_seq = (
        orig_seq[:interval_start] +
        str(NON_CODING_TRANSITION_TABLE[orig_seq[interval_start]]) +
        orig_seq[interval_start + 1:])
    assert len(orig_seq) == len(updated_seq)
    mutable_genome_record.seq = updated_seq
    return mutable_genome_record


def _remove_site_in_coding_feature(
        refactor_context,
        mutable_genome_record,
        site_occur_obj,
        feature):
    """Removes a restriction site that falls in a region annotated as
    coding.

    The strategy is to muddle up all affected codons within the feature
    so as to reduce the chance of "snap-back" over generations.

    Args:
        refactor_context: The RefactorContext.
        mutable_genome_record: The SeqRecord object representing the genome.
        site_occur_obj: Object containing data about a specific occurrence of a
            restriction enzyme site in the genome.
        feature: The feature that is overlapped by the restriction site.

    Returns:
        Object with keys:
            * is_success: Whether remove succeeded.
            * updated_genome_record: The updated genome_record if successful.
            * exception_string: Message describing failure.
    """
    interval = site_occur_obj['interval']

    # Figure out the specific codons that need to be changed.
    affected_codon_indeces = get_region_codon_indeces_in_feature(
            feature, interval)
    avoid_codons_in_positions = {}
    feature_seq = str(feature.extract(mutable_genome_record.seq))
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


def _remove_site_with_partial_feature_overlap(
        mutable_genome_record, site_occur_obj, feature):
    """DEPRECATED.

    Removes a restriction site that partially overlaps a feature.

    The strategy is to identify the end of the restriction site that is
    not inside of the feature and mutate the first base on that end
    according to NON_CODING_TRANSITION_TABLE.
    """
    orig_seq = mutable_genome_record.seq

    interval = site_occur_obj['interval']
    interval_start = interval[0]
    interval_end = interval[1]

    if feature.location.start <= interval_start:
        # >>>>>
        #    (....)
        assert interval_end > feature.location.end
        updated_seq = (
            orig_seq[:interval_end - 1] +
            str(NON_CODING_TRANSITION_TABLE[orig_seq[interval_end - 1]]) +
            orig_seq[interval_end:])
    else:
        #    >>>>>
        # (....)
        updated_seq = (
            orig_seq[:interval_start] +
            str(NON_CODING_TRANSITION_TABLE[orig_seq[interval_start]]) +
            orig_seq[interval_start + 1:])

    assert len(orig_seq) == len(updated_seq)
    mutable_genome_record.seq = updated_seq
    return mutable_genome_record


if __name__ == '__main__':
    import cProfile
    from datetime import datetime
    import os
    import pickle

    from Bio import SeqIO

    from biopython_util import get_genome_record
    from refactor_config import OUTPUT_DIR
    from refactor_config import RESTRICTION_ENZYME_SITES_TO_REMOVE
    from refactor_context import RefactorContext

    genome_record = get_genome_record(
            '../data/2013_03_06_20_16_04_mds42_refactored.gbk')

    START_BOUND = 100000
    END_BOUND = 150000

    for enzyme in RESTRICTION_ENZYME_SITES_TO_REMOVE:
        runs = find_restriction_site_occurrences(
                genome_record, enzyme, start_bound=START_BOUND,
                end_bound=END_BOUND)
        print enzyme
        print runs

    # refactor_context = RefactorContext(genome_record)

    # restriction_site_list = RESTRICTION_ENZYME_SITES_TO_REMOVE

    # print 'Removing sites...'
    # tmp_file_prefix = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    # cprofile_output_dest = os.path.join(
    #         OUTPUT_DIR, tmp_file_prefix + '_cprofile.out')
    # cmd = ('result = remove_restriction_sites('
    #         'refactor_context, restriction_site_list,'
    #         'start_bound=START_BOUND, end_bound=END_BOUND)')
    # cProfile.run(cmd, cprofile_output_dest)

    # updated_genome_record = result['updated_genome_record']
    # flagged = result['flagged']

    # # Write the resulting genome.
    # genome_output_file = (
    #         '../data/2013_03_06_20_16_04_mds42_refactored_res_sites_removed.gbk'
    # )
    # with open(genome_output_file, 'w') as output_fh:
    #     SeqIO.write(updated_genome_record, output_fh, 'genbank')

    # # Write the flagged data.
    # metadata_output_file = (
    #         '../data/2013_03_06_20_16_04_mds42_refactored_res_sites_removed.flagged'
    # )
    # with open(metadata_output_file, 'w') as metadata_output_fh:
    #     pickle.dump(flagged, metadata_output_fh)
