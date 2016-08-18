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
Methods for performing various analysis related to the project.
"""

import csv
import json
import os
import pickle

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
import pandas as pd

# Try to use GETK for future development.
from getk.biopython_util import get_feature_by_gene_name
from getk.biopython_util import get_feature_gene_name
from getk.biopython_util import get_feature_qualifier_value
from getk.biopython_util import get_global_position
from getk import codon_usage_memex
from getk.codon_usage_memex import get_ecoli_codon_usage_memex


# Legacy imports

import biopython_util
from biopython_util import get_genome_record
from biopython_util import get_feature_by_id
from biopython_util import get_feature_gene
from biopython_util import does_interval_overlap_feature
from conflicting_pair_finder import find_all_overlaps
from homology import calc_homology
from homology import find_features_to_check_for_homology
from refactor_config import CACHE_DIR
from refactor_config import CODONS_TO_REMOVE
from refactor_config import RESTRICTION_ENZYME_SITES_TO_REMOVE
from paths import DATA_DIR
from paths import GENOMES_DIR
from per_base_conservation import get_mds42_per_base_conservation
from post_processing.gc_content_fixer import find_gc_content_extremes
from post_processing.homopolymers import find_homopolymer_runs
from post_processing.restriction_sites import find_restriction_site_occurrences
from tool_wrappers.rbs_calc_util import calc_rbs_score_for_feature
from tool_wrappers.rbs_calc_util import get_mds42_rbs_strength_profile
from tool_wrappers.rbs_calc_util import RBS_CALC_BUFFER


# These were manually modified to preserve SECIS sites.
MANUAL_SECIS_FIXES = set([
    'fdhF',
    'fdnG',
    'fdoG',
])

MANUAL_AGR_FIXES = set([
    'secE',
])

FORBIDDEN_CODONS = set([
    'TTA',
    'TTG',
    'AGA',
    'AGG',
    'AGC',
    'AGT',
    'TAG'
])


def determine_feature_frame(feature):
    """Determine the gene reading frame, relative to the first base in the
    genome.

    Args:
        feature: Bio.SeqFeature object.

    Returns:
        An element in {-3, -2, -1, 1, 2, 3}. For example, a result of 1 would
        mean that the gene would be transcribed in the same frame as if the
        whole genome were transcribed from the first base.
    """
    frame = feature.location.start % 3 + 1
    if feature.strand == -1:
        frame *= -1
    return frame


def determine_codon_usage_table(
        genome_record,
        codon_usage_source_file='../data/ecoli-codon-usage.txt',
        ignore_feature_ids=[]):
    """Look at all the CDS features in a genome and determine the codon usage
    distribution.

    Args:
        genome_record: A SeqRecord with features representing genomes
        codon_usage_source_file: Source file used to determine amino acid
            to codon mapping.

    Returns:
        A dictionary with keys being codons and values being dictionaries
        containing data about the codon, keyed by the data label, including:
            * usage
            * count
            * amino_acid
            * amino_acid_count
    """
    # Next create a table from amino acid to codon. We use the source
    # file to get the mapping from amino acid to codon, but we will still
    # handle calculating the distribution ourselves below.
    data_as_strings = []
    with open(codon_usage_source_file) as fh:
        for line in fh:
            splitline = line.split()
            for i in [0, 5, 10, 15]:
                data_as_strings.append(splitline[i:i+5])

    # Build a dictionary from amino acid to list of codons.
    print 'Building the dictionaries from source file...'
    aa_to_codon_dict = {}
    codon_usage_dict = {}
    for string in data_as_strings:
        amino_acid = string[1]
        dna_codon = string[0].replace('U','T')
        if not amino_acid in aa_to_codon_dict:
            aa_to_codon_dict[amino_acid] = {}
        aa_to_codon_dict[amino_acid][dna_codon] = {}
        if not dna_codon in codon_usage_dict:
            codon_usage_dict[dna_codon] = {
                    'amino_acid': amino_acid,
                    'count': 0,
                    'usage': 0.0,
            }

    # We only use CDS features for this stat.
    coding_features = filter(lambda feature: feature.type == 'CDS',
            genome_record.features)

    # Ignore removed.
    # NOTE: Slightly awkward due to id format having been changed
    # since this data run.
    len_before_ignore = len(coding_features)
    prefix_ignore_ids = map(lambda x: x[:15], ignore_feature_ids)
    coding_features = filter(
            lambda feature: not feature.id[:15] in prefix_ignore_ids,
            coding_features)
    len_after_ignore = len(coding_features)
    assert len_after_ignore == len_before_ignore - len(ignore_feature_ids)

    genome_record_seq_str = str(genome_record.seq)

    # First count up the usage of all codons.
    print 'Calculating running count...'
    codon_count = {}
    for f_idx, feature in enumerate(coding_features):
        # DEBUG
        # if f_idx % 100 == 0:
        #     print 'Counting %d of %d' % (f_idx, len_after_ignore)

        seq = feature.extract(genome_record_seq_str)
        assert len(seq) % 3 == 0
        for codon_index in range(0, len(seq), 3):
            codon = str(seq[codon_index:codon_index + 3])
            if not codon in codon_count:
                codon_count[codon] = 1
            else:
                codon_count[codon] += 1

    # Now calculate the distributions and create the object to return.
    print 'Calculating distributions...'
    for amino_acid, codon_dict in aa_to_codon_dict.iteritems():
        # Caculate the denominator.
        total_codons_for_amino_acid = 0
        for codon in codon_dict.keys():
            total_codons_for_amino_acid += codon_count.get(codon, 0)

        # Caculate the ratio and add it to the result dictionary.
        for codon in codon_dict.keys():
            codon_usage_dict[codon]['count'] = codon_count.get(codon, 0)
            codon_usage_dict[codon]['amino_acid_count'] = (
                    total_codons_for_amino_acid)
            codon_usage_dict[codon]['usage'] = (
                    float(codon_count.get(codon, 0)) /
                            total_codons_for_amino_acid
            )

    print 'Done..'

    return codon_usage_dict


def find_non_terminating_stop_codons_in_feature(genome_record):
    STOP_CODONS = set(['TAA', 'TGA', 'TAG'])

    found = {}

    # We only use CDS features for this stat.
    coding_features = filter(lambda feature: feature.type == 'CDS',
            genome_record.features)

    for feature in coding_features:
        seq = feature.extract(genome_record).seq
        for codon_start_pos in range(0, len(seq) - 3, 3): # Ignore last one.
            codon = str(seq[codon_start_pos:codon_start_pos+ 3])
            if codon in STOP_CODONS:
                if not feature.id in found:
                    found[feature.id] = []
                found[feature.id].append({
                    'codon_start_pos': codon_start_pos,
                    'mds42_position':
                            feature.location.start + codon_start_pos + 1,
                })

    return found


def determine_hexamer_usage_dict(genome_record, dicodons_only=False):
    """Step through each coding feature, one base at a time, and keep track
    of the number of hexamers. Finally, report the proportion of that hexamer
    as a fraction of all hexamers observed.
    """
    # We only use CDS features for this stat.
    coding_features = filter(lambda feature: feature.type == 'CDS',
            genome_record.features)

    # Iterate in codon-step sizes if only want dicodons.
    if dicodons_only:
        iter_step = 3
    else:
        iter_step = 1

    hexamer_count  = {}
    for feature in coding_features:
        seq = feature.extract(genome_record).seq
        # Iterate within feature only.
        for base_index in range(0, len(seq) - 5, iter_step):
            hexamer = str(seq[base_index:base_index + 6])
            if not hexamer in hexamer_count:
                hexamer_count[hexamer] = 1
            else:
                hexamer_count[hexamer] += 1

    total_hexamers = sum(hexamer_count.values())

    # Calculate distributions.
    hexamer_dist = {}
    for hexamer, count in hexamer_count.iteritems():
        hexamer_dist[hexamer] = {
                'count': count,
                'usage': float(count) / total_hexamers,
        }

    return hexamer_dist


def determine_dicodon_usage_dict(genome_record):
    """Calculate dicodon usage.
    """
    return determine_hexamer_usage_dict(genome_record, dicodons_only=True)


def get_failed_genes_from_metadata(meta_file_location):
    failed = []
    with open(meta_file_location) as fh:
        meta = pickle.load(fh)
    for feature_id, data in meta.iteritems():
        if not data['is_success']:
            failed.append((feature_id, data))
    return failed


def do_genbank_gene_features_contain_more_data_than_respective_coding_element():
    """One-time use function used to figure out whether genbank 'gene' feature
    annotations contain any data that their corresponding coding element doesn't.
    I'm wondering about this to make sure that it's safe to ignore features while
    recoding, and then just add them back in as a subset of their child object.
    """
    record = biopython_util.get_genome_record(
            '../data/mds42_full.gbk', features_to_ignore=['source'])


    # We've previously established that there are two features for each
    # locus_tag: One is of type 'gene' and the other is of another type such
    # such as 'CDS', 'misc_feature', etc.
    map_locus_tag_to_annotation_pair = {}
    for feature in record.features:
        locus_tag = feature.qualifiers['locus_tag'][0]
        pair = map_locus_tag_to_annotation_pair.get(locus_tag, [None, None])
        if feature.type == 'gene':
            assert not pair[0]
            pair[0] = feature
        else:
            assert not pair[1]
            pair[1] = feature
        map_locus_tag_to_annotation_pair[locus_tag] = pair

    # Make sure above step worked.
    for locus_tag, pair in map_locus_tag_to_annotation_pair.iteritems():
        assert 2 == len(pair) and pair[0] and pair[1]
    print 'Asserted that locus_tags come in pairs.'

    # Now check that feature qualifiers of gene is subset of feature qualifiers
    # of corresponding coding sequence in the genbank file.
    for pair in map_locus_tag_to_annotation_pair.values()[:10]:
        assert set(pair[0].qualifiers.keys()).issubset(pair[1].qualifiers.keys())
    print 'Successfully confirmed that all gene feature qualifiers are subsets.'


def add_rbs_scores_to_feature_qualifiers(refactored_genome_record):
    """Scan through the features in a genome record and add original and new
    RBS scores to the Biopython feature.qualifiers dictionary (which also
    shows up as the metadata for each feature in the genbank file).

    NOTE: Works for MDS42 default only right now.
    """
    genome_forward_seq = refactored_genome_record.seq

    features_to_profile = filter(
            lambda feature: feature.type in ['CDS', 'misc_feature'],
            refactored_genome_record.features
    )

    orig_rbs_profile = get_mds42_rbs_strength_profile()

    num_features = len(features_to_profile)
    for feature_index in range(num_features):
        feature = features_to_profile[feature_index]

        print 'Calculating feature %d of num_features %d' % (
                    feature_index + 1, num_features)
        print '...Feature id: %s' % feature.id

        # Get the original expression from the pre-calculated table.
        orig_expression = orig_rbs_profile[feature.id]

        # Calculate the new expression score.
        refactored_expression = calc_rbs_score_for_feature(
                feature, genome_forward_seq, RBS_CALC_BUFFER)['expression']

        # Cast to ints for saving.
        if isinstance(orig_expression, float):
            orig_expression = int(orig_expression)
        if isinstance(refactored_expression, float):
            refactored_expression = int(refactored_expression)

        feature.qualifiers['rbs_expr_orig'] = orig_expression
        feature.qualifiers['rbs_expr_refactored'] = refactored_expression
        print '...original expression:', feature.qualifiers['rbs_expr_orig']
        print '...refactored expression:', feature.qualifiers['rbs_expr_refactored']


    with open(refactored_genbank_source, 'w') as output_fh:
        SeqIO.write(refactored_genome_record, output_fh, 'genbank')


def calc_homology_for_copy_fixes_in_genome_record(
        refactored_genbank_source, outfile):
    """For features in a refactored genome that are created as a result of
    copying another part, this function calculates the homology between
    the copy and the original part.

    This method creates a report, using the logic for measuring homology
    in homology.py. That module also contains logic for fixing homology.

    The results are written to a file.
    """
    refactored_genome_record = biopython_util.get_genome_record(
            refactored_genbank_source)

    homologous_pair_obj_list = find_features_to_check_for_homology(
            refactored_genome_record)

    # Calculate scores for all of them.
    updated_homologous_pairs = []
    for pair_obj in homologous_pair_obj_list:
        updated_pair_obj = pair_obj.copy()
        if pair_obj['source_seq'] == 'UNKNOWN' or pair_obj['copy_seq'] == 'UNKNOWN':
            updated_pair_obj['homology'] = 'UNKNOWN'
        else:
            updated_pair_obj['homology'] = calc_homology(pair_obj)
        updated_homologous_pairs.append(updated_pair_obj)

    # Write the output to a .csv file.
    FIELD_NAMES = [
        'type',
        'source_id',
        'source_start',
        'homology',
        'source_seq_first_codon_index',
        'source_seq',
        'copy_id',
        'copy_seq',
    ]

    # Write the result to file.
    with open(outfile, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, FIELD_NAMES)
        writer.writeheader()

        for pair_obj in updated_homologous_pairs:
            writer.writerow(pair_obj)


def analyze_restriction_enzyme_sites():
    """Analyze where restriction enzyme sites fall.
    """
    genome_record = get_genome_record('../data/mds42_full.gbk')

    # Write the output to a .csv file.
    FIELD_NAMES = [
        'site',
        'interval',
        'strand',
        'overlaps_feature_id',
        'overlaps_feature_type',
        'overlaps_feature_strand',
    ]

    outfile = '../data/mds42_restriction_enzyme_analysis.csv'
    with open(outfile, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, FIELD_NAMES)
        writer.writeheader()

        for site in RESTRICTION_ENZYME_SITES_TO_REMOVE:
            site_occurrences = find_restriction_site_occurrences(
                    genome_record, site)
            for site_occurrence in site_occurrences:
                interval = site_occurrence['interval']
                strand = site_occurrence['strand']
                found_overlap = False
                last_feature = None
                for feature in genome_record.features:
                    overlap = does_interval_overlap_feature(interval, feature)
                    if overlap:
                        # if found_overlap:
                        #     assert False, "site_obj:\n%s\n\nlast:\n%s\n\nthis:\n%s" % (str(site_occurrence), str(last_feature), str(feature))
                        writer.writerow({
                            'site': str(site),
                            'interval': interval,
                            'strand': strand,
                            'overlaps_feature_id': feature.id,
                            'overlaps_feature_type': feature.type,
                            'overlaps_feature_strand': feature.strand
                        })
                        found_overlap = True
                        break
                    last_feature = feature
                if not found_overlap:
                    writer.writerow({
                        'site': str(site),
                        'interval': interval,
                        'strand': strand,
                        'overlaps_feature_id': 'NONE',
                        'overlaps_feature_type': 'NONE',
                        'overlaps_feature_strand': 'NONE'
                    })


def find_rna_genes(input_genbank_source, out_csv_file):
    """Identifies the list of RNA genes in a Genome Record and writes the
    result to a .csv file.
    """
    RNA_TYPES = set([
        'ncRNA',
        'rRNA',
        'tRNA',
        'tmRNA'
    ])

    FIELD_NAMES = [
        'gene',
        'rna_type'
    ]

    genome_record = get_genome_record(
            input_genbank_source, assert_unique_ids=False)

    with open(out_csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, FIELD_NAMES)
        writer.writeheader()

        for feature in genome_record.features:
            if feature.type in RNA_TYPES:
                writer.writerow({
                    'gene': feature.qualifiers['gene'][0],
                    'rna_type': feature.type
                })


def generate_gene_start_end_table(input_genbank_source, out_csv_file,
        pythonic_indeces=True):
    """Generates a table with gene name, start, end, polarity."""
    FIELD_NAMES = [
        'gene',
        'gene_start',
        'gene_end',
        'gene_strand',
    ]

    genome_record = get_genome_record(
            input_genbank_source,
            only_get_features=['gene'],
            assert_unique_ids=False)

    with open(out_csv_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, FIELD_NAMES)
        writer.writeheader()

        for feature in genome_record.features:
            if 'gene' in feature.qualifiers:
                start = feature.location.start
                end = feature.location.end
                if pythonic_indeces:
                    start += 1
                    end += 1
                writer.writerow({
                    'gene': feature.qualifiers['gene'][0],
                    'gene_start': start,
                    'gene_end': end,
                    'gene_strand': feature.strand,
                })


def generate_recoding_stats(original_genome_record, recoded_genome_record,
        outfile):
    """Generate stats related to recoding.
    """
    # We only use CDS features for this stat.
    original_coding_features = filter(
            lambda feature: feature.type == 'CDS',
            original_genome_record.features)

    total_coding_features = len(original_coding_features)
    total_codons = 0
    num_codons_recoded = 0
    percent_recoded_list = []
    for idx, feature in enumerate(original_coding_features):
        print 'Analyzing feature %s,  %d of %d' % (
                feature.id, idx + 1, total_coding_features)
        recoded_feature = get_feature_by_id(recoded_genome_record, feature.id)
        if recoded_feature is None:
            continue
        orig_seq = feature.extract(original_genome_record).seq
        recoded_seq = recoded_feature.extract(recoded_genome_record).seq
        codons_recoded_for_feature = 0
        feature_len = len(orig_seq)
        for codon_index in range(0, len(orig_seq), 3):
            total_codons += 1
            original_codon = str(orig_seq[codon_index:codon_index + 3])
            recoded_codon = str(recoded_seq[codon_index:codon_index + 3])
            if original_codon != recoded_codon:
                codons_recoded_for_feature += 1
        feature_codons = feature_len / 3
        feature_percent_recoded = (float(codons_recoded_for_feature) /
                feature_codons)
        percent_recoded_list.append(feature_percent_recoded)

        num_codons_recoded += codons_recoded_for_feature

    percent_recoded = float(num_codons_recoded) / total_codons
    avg_percent_recoded = (float(sum(percent_recoded_list)) /
            len(percent_recoded_list))

    data = {
        'total_CDS_features': total_coding_features,
        'total_codons': total_codons,
        'num_codons_recoded': num_codons_recoded,
        'proportion_recoded': percent_recoded,
        'avg_proportion_recoded': avg_percent_recoded
    }

    with open(outfile, 'w') as fh:
        fh.write(json.dumps(data, indent=4, separators=(',', ': ')))


def count_forbidden_codons(genome_record):
    count = 0
    total_codons = 0
    forbidden_codons = set(CODONS_TO_REMOVE)

    original_coding_features = filter(
            lambda feature: feature.type == 'CDS',
            genome_record.features)

    total_coding_features = len(original_coding_features)
    for idx, feature in enumerate(original_coding_features):
        print 'Analyzing feature %s,  %d of %d' % (
                feature.id, idx + 1, total_coding_features)
        orig_seq = feature.extract(genome_record).seq
        for codon_index in range(0, len(orig_seq), 3):
            total_codons += 1
            original_codon = str(orig_seq[codon_index:codon_index + 3])

            if original_codon in forbidden_codons:
                count += 1

    print count, total_codons


def analyze_conservation_data_utility_in_trouble_regions(genome_record,
        outfile):
    """Analyzes the genome record to see how useful the conservation data is
    for the given regions.
    """
    print 'Getting trouble regions...'
    trouble_regions = identify_trouble_regions(genome_record)

    # Get each individual position.
    trouble_positions = set()
    for interval_tuple in trouble_regions:
        trouble_positions |= set(range(*interval_tuple))

    # Ignore coding regions.
    ignore_coding_regions = [(feature.location.start, feature.location.end)
            for feature in genome_record.features if feature.type == 'CDS']
    for interval_tuple in ignore_coding_regions:
        trouble_positions -= set(range(*interval_tuple))

    aggregate_result = get_conservation_aggregate_for_positions(
            trouble_positions, genome_record)
    at_least_one_alternate = aggregate_result['at_least_one_alternate']
    greater_than_one_alternate = aggregate_result['greater_than_one_alternate']

    data = {
        'total_trouble_positions': len(trouble_positions),
        'at_least_one_alternates': at_least_one_alternate,
        'greater_than_one_alternate': greater_than_one_alternate
    }
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(data, indent=4, separators=(',', ': ')))


def get_conservation_aggregate_for_positions(positions, genome_record,
        mds42_conservation_dict=None):
    """Performs aggregate count of conservation over positions.

    Returns:
        A dict with keys:
            * at_least_one_alternate
            * greater_than_one_alternate
    """
    if not mds42_conservation_dict:
        print 'Getting conservation dict...'
        mds42_conservation_dict = get_mds42_per_base_conservation()

    print 'Calculating aggregate...'
    at_least_one_alternate = 0
    greater_than_one_alternate = 0
    BASES = set(['A', 'C', 'G', 'T'])
    for position in positions:
        at_least_one = False
        greater_than_one = False
        conservation_data = mds42_conservation_dict[position]
        num_genomes = conservation_data['num_genomes']
        orig_base = genome_record.seq[position]
        remaining_bases = BASES - set(orig_base)
        for base in remaining_bases:
            occurrences = int(conservation_data[base] * num_genomes)
            if occurrences > 0:
                at_least_one = True
            if occurrences > 1:
                greater_than_one = True
        if at_least_one:
            at_least_one_alternate += 1
        if greater_than_one:
            greater_than_one_alternate += 1

    return {
        'at_least_one_alternate': at_least_one_alternate,
        'greater_than_one_alternate': greater_than_one_alternate
    }


def identify_trouble_regions(genome_record, force_recalculate=False):
    """Identify trouble regions in the genome.

    Our current definition of trouble regions is:
        (1) Occur outside of a gene AND any of:
            * GC content extremes
            * Homopolymer runs
            * Restriction site

    Returns:
        List of two-tuples representing pythonic intervals in genome_record.
    """
    TROUBLE_REGIONS_CACHE = os.path.join(CACHE_DIR, 'trouble_regions.cache')

    trouble_regions = []

    if not force_recalculate:
        try:
            with open(TROUBLE_REGIONS_CACHE) as fh:
                return pickle.load(fh)
        except:
            # Cache doesn't exist. Add more precise error catching if desired.
            pass

    # GC content
    print '...Analyzing GC content...'
    trouble_regions += [res['interval'] for res in
            find_gc_content_extremes(genome_record)]

    # Homopolymer runs
    print '...Analyzing homopolymer runs...'
    trouble_regions += [res['interval'] for res in
            find_homopolymer_runs(genome_record)]

    # Restriction sites.
    print '...Analyzing restriction sites...'
    for res_site in RESTRICTION_ENZYME_SITES_TO_REMOVE:
        trouble_regions += [res['interval'] for res in
                find_restriction_site_occurrences(genome_record, res_site)]

    print '...Caching the result...'
    with open(TROUBLE_REGIONS_CACHE, 'w') as fh:
        pickle.dump(trouble_regions, fh)
    print '...Done caching.'

    return trouble_regions


def analyze_conservation_in_genes(genome_record, outfile):
    """Analyze conservation data in genes to get a sense of how much
    conservation there is.
    """
    coding_regions = set([(feature.location.start, feature.location.end)
            for feature in genome_record.features if feature.type == 'CDS'])

    all_cds_positions = set()
    for interval_tuple in coding_regions:
        all_cds_positions |= set(range(*interval_tuple))

    aggregate_result = get_conservation_aggregate_for_positions(
            all_cds_positions, genome_record)
    at_least_one_alternate = aggregate_result['at_least_one_alternate']
    greater_than_one_alternate = aggregate_result['greater_than_one_alternate']

    data = {
        'total_cds_positions': len(all_cds_positions),
        'cds_positions_with_at_least_one_alternate': at_least_one_alternate,
        'cds_positions_with_greater_than_one_alternate':
                greater_than_one_alternate
    }
    with open(outfile, 'w') as fh:
        fh.write(json.dumps(data, indent=4, separators=(',', ': ')))


def analyze_variation(original_record, recoded_record):
    """Analyzes the variation between the two records.

    The start and recoded records must have the same feature ids.

    Returns:
        Dictionary with the following keys:
            * codon_similarity: Proportion of codons unchanged.
            * based_similarity: Proportion of bases unchanged.
    """
    total_codons = 0
    reassigned_codons = 0
    total_bases = 0
    reassigned_bases = 0

    original_record_coding_features = [feature for feature in
            original_record.features if feature.type == 'CDS']

    for orig_feature in original_record_coding_features:
        recoded_feature = get_feature_by_id(recoded_record, orig_feature.id)

        orig_seq = orig_feature.extract(original_record.seq)
        recoded_seq = recoded_feature.extract(recoded_record.seq)

        if not len(orig_seq) == len(recoded_seq):
            print '>>>> Omitting ' + str(orig_feature)
            continue

        total_bases += len(orig_feature)
        total_codons += len(orig_feature) / 3

        for codon_index in range(0, len(orig_seq), 3):
            orig_codon = str(orig_seq[codon_index:codon_index + 3])
            recoded_codon = str(recoded_seq[codon_index:codon_index + 3])
            if orig_codon == recoded_codon:
                continue
            else:
                reassigned_codons += 1
                reassigned_bases += _num_differing_bases(orig_codon,
                        recoded_codon)

    return {
        'total_codons': total_codons,
        'reassigned_codons': reassigned_codons,
        'codon_similarity': float(total_codons - reassigned_codons) / total_codons,
        'total_bases': total_bases,
        'reassigned_bases': reassigned_bases,
        'base_similarity': float(total_bases - reassigned_bases) / total_bases
    }


def _num_differing_bases(seq1, seq2):
    """Return the number of different bases between the two sequences.
    """
    assert len(seq1) == len(seq2)
    num_different = 0
    for idx, seq1_base in enumerate(seq1):
        seq2_base = seq2[idx]
        if seq1_base != seq2_base:
            num_different += 1
    return num_different


def analyze_codon_motif(motif_set, original_record, recoded_record, report_file):
    """Find all AGY positions that are are within 30 bases upstream of another
    feature.

    Args:
        motif_set: Set of codons to look for.
    """
    AGY_list = []

    original_seq = str(original_record.seq).upper()
    recoded_seq = str(recoded_record.seq).upper()

    # Identify overlapping features (this includes those close enough for RBS)
    # to possibly conflict. Use cache.
    conflicting_pairs = find_all_overlaps(
        None, # genome_record,
        None, # forbidden_codons,
        cache=True)

    original_rbs_profiles = get_mds42_rbs_strength_profile()

    for pair in conflicting_pairs:
        if (pair['upstream_feature'].type != 'CDS' or
                pair['downstream_feature'].type != 'CDS'):
            continue
        # Ignore prfB.
        if (get_feature_gene(pair['upstream_feature']) == 'prfB' or
                get_feature_gene(pair['downstream_feature']) == 'prfB'):
            continue
        original_upstream_feature = get_feature_by_id(
                original_record, pair['upstream_feature'].id)
        original_downstream_feature = get_feature_by_id(
                original_record, pair['downstream_feature'].id)
        recoded_upstream_feature = get_feature_by_id(
                recoded_record, pair['upstream_feature'].id)
        recoded_downstream_feature = get_feature_by_id(
                recoded_record, pair['downstream_feature'].id)

        if (original_upstream_feature.strand == 1 and
                original_downstream_feature.strand == 1):
            roi_start = original_downstream_feature.location.start - 20
            for codon_index in range(0, len(original_upstream_feature), 3):
                pos = original_upstream_feature.location.start + codon_index
                if pos < roi_start:
                    continue
                else:
                    codon = original_seq[pos:pos +3]
                    if codon in motif_set:
                        recoded_pos = (recoded_upstream_feature.location.start +
                                codon_index)
                        recoded_codon = recoded_seq[recoded_pos:recoded_pos +3]

                        # Get before/after RBS.
                        orig_rbs_expression = original_rbs_profiles[
                                original_downstream_feature.id]
                        recoded_rbs_expression = calc_rbs_score_for_feature(
                                recoded_downstream_feature, recoded_seq)[
                                        'expression']
                        delta_expression = (recoded_rbs_expression -
                                orig_rbs_expression)

                        AGY_list.append({
                            'pos': pos,
                            'recoded_pos': recoded_pos,
                            'ref': codon,
                            'alt': recoded_codon,
                            'orig_rbs_expression': orig_rbs_expression,
                            'recoded_rbs_expression': recoded_rbs_expression,
                            'delta_expression': delta_expression,
                            'strand': 1,
                        })

        # elif (original_upstream_feature.strand == -1 and
        #         original_downstream_feature.strand == -1):
        #     original_downstream_feature_seq = (
        #             original_downstream_feature.extract(original_seq))
        #     recoded_downstream_feature_seq = (
        #             recoded_downstream_feature.extract(recoded_seq))
        #     roi_start = original_upstream_feature.location.end + 20
        #     for codon_index in range(0, len(original_downstream_feature), 3):
        #         pos = original_downstream_feature.location.end - codon_index
        #         if pos > roi_start:
        #             continue
        #         else:
        #             codon = original_downstream_feature_seq[
        #                     codon_index:codon_index + 3]
        #             if codon in AGY_codons:
        #                 recoded_pos = (recoded_downstream_feature.location.start -
        #                         codon_index)
        #                 recoded_codon = recoded_downstream_feature_seq[
        #                         codon_index:codon_index + 3]

        #                 # Get before/after RBS.
        #                 orig_rbs_expression = original_rbs_profiles[
        #                         original_upstream_feature.id]
        #                 recoded_rbs_expression = calc_rbs_score_for_feature(
        #                         recoded_upstream_feature, recoded_seq)[
        #                                 'expression']
        #                 delta_expression = (recoded_rbs_expression -
        #                         orig_rbs_expression)

        #                 AGY_list.append({
        #                     'pos': pos,
        #                     'recoded_pos': recoded_pos,
        #                     'ref': codon,
        #                     'alt': recoded_codon,
        #                     'orig_rbs_expression': orig_rbs_expression,
        #                     'recoded_rbs_expression': recoded_rbs_expression,
        #                     'delta_expression': delta_expression,
        #                     'strand': -1
        #                 })

    print 'Writing report.'
    with open(report_file, 'w') as fh:
        FIELD_NAMES = [
            'pos',
            'recoded_pos',
            'ref',
            'alt',
            'orig_rbs_expression',
            'recoded_rbs_expression',
            'delta_expression',
            'strand',
        ]

        writer = csv.DictWriter(fh, FIELD_NAMES)
        writer.writeheader()
        for AGY_obj in AGY_list:
            writer.writerow(AGY_obj)


def build_gene_to_segment_map(rc3_seq_record):
    """Buids map from gene to corresponding segment in annotation.
    """
    segment_features = [f for f in rc3_seq_record.features
            if f.type == 'synthesis_seg']
    segment_features = sorted(segment_features, key=lambda f: f.location.start)

    cds_features = [f for f in rc3_seq_record.features if f.type == 'CDS']
    cds_features = sorted(cds_features, key=lambda f: f.location.start)
    num_cds_features = len(cds_features)

    cds_features_idx = 0
    map_gene_to_segment = {}
    for seg_f in segment_features:
        seg_label = get_feature_qualifier_value(seg_f, 'label')
        cds_f = cds_features[cds_features_idx]
        while cds_f.location.start < seg_f.location.end:
            map_gene_to_segment[get_feature_gene_name(cds_f)] = seg_label
            cds_features_idx += 1
            if cds_features_idx >= num_cds_features:
                break
            cds_f = cds_features[cds_features_idx]
    assert len(cds_features) == len(map_gene_to_segment), (
            "Expected: %d. Mapped: %d" %
                    (len(cds_features), len(map_gene_to_segment)))

    return map_gene_to_segment


def generate_list_of_recoded_codons(
        orig_seq_record, recoded_seq_record, output_report,
        annotations_fixed_recoded_seq_record):
    """Analyze RC3 codon scores.

    Generates table with the following common keys:
        * segment
        * gene
        * strand
        * codon_index
        * position_mds42
        * position_rc3
        * codon_mds42
        * codon_rc3

    NOTE: Copied from getk_validation/scripts/analyze_rc3.py#analyze_all_codons
    on 2016-01-13. I had previously used that code to generate this, but that
    code has since changed so just reworking it here so we have it.
    """
    orig_seq = str(orig_seq_record.seq)
    recoded_seq = str(recoded_seq_record.seq)

    codon_usage_memex_obj = get_ecoli_codon_usage_memex(randomized=True)

    START_CODON_SET = set(codon_usage_memex.CodonUsageMemex.start_codons)

    data_obj_list = []

    # For segment annotation.
    map_gene_to_segment = build_gene_to_segment_map(recoded_seq_record)

    orig_cds_features = [f for f in orig_seq_record.features if f.type == 'CDS']

    # Manually include prfB gene. No CDS in original MDS42 record because of
    # frameshift.
    orig_cds_features.append(
            get_feature_by_gene_name('prfB', 'gene', orig_seq_record))
    # orig_cds_features = [get_feature_by_gene_name('prfB', 'gene', orig_seq_record)]
    orig_cds_features = sorted(
            orig_cds_features, key=lambda f: f.location.start)

    total_coding_features = len(orig_cds_features)

    # Convenience map to make gene lookups below constant time.
    map_gene_name_to_feature_in_recoded = {}
    for f in recoded_seq_record.features:
        if not f.type == 'CDS':
            continue
        gene_name = get_feature_gene_name(f)
        assert not gene_name in map_gene_name_to_feature_in_recoded, gene_name
        map_gene_name_to_feature_in_recoded[gene_name] = f

    # Annotated forbidden codon positions.
    forbidden_codon_positions = set()
    for f in annotations_fixed_recoded_seq_record.features:
        if f.type == 'recoded_codon':
            if f.strand == 1:
                forbidden_codon_positions.add(int(f.location.start))
            else:
                forbidden_codon_positions.add(int(f.location.start) + 2)
    assert len(forbidden_codon_positions) == 62214, len(forbidden_codon_positions)

    # Track codon positions we've seen and compare to forbiddenCodon
    # annotations.
    observed_forbidden_codon_positions = set()

    for idx in range(0, total_coding_features):
        feature = orig_cds_features[idx]

        # DEBUG: Verbose output.
        # if idx % 100 == 0:
        #     print 'Analyzing %d-strand feature %s,  %d of %d' % (
        #             feature.strand, get_feature_gene_name(feature),
        #             idx + 1, total_coding_features)

        gene_name = get_feature_gene_name(feature)
        assert gene_name is not None

        recoded_f = map_gene_name_to_feature_in_recoded.get(gene_name)
        if recoded_f is None:
            print 'Gene %s (%d) not found in recoded genome. Skipping.' % (
                    gene_name, feature.location.start)
            continue

        if gene_name == 'prfB':
            assert len(feature) == len(recoded_f) + 1
        else:
            assert len(feature) == len(recoded_f)

        orig_feature_seq = feature.extract(orig_seq)
        recoded_feature_seq = recoded_f.extract(recoded_seq)

        # NOTE: Skip first codon. E.g. sometimes this is TTG.
        codon_index = 3
        while codon_index < len(feature):
            data_obj = {}

            # Manually frameshift.
            maybe_prfB_hack_codon_index = codon_index
            if gene_name == 'prfB' and codon_index > 72:
                maybe_prfB_hack_codon_index += 1

            orig_codon = orig_feature_seq[
                    maybe_prfB_hack_codon_index:maybe_prfB_hack_codon_index + 3]
            if not orig_codon in FORBIDDEN_CODONS:
                codon_index += 3
                continue

            recoded_codon = recoded_feature_seq[codon_index:codon_index + 3]
            try:
                assert codon_usage_memex_obj.are_synonymous(
                        orig_codon, recoded_codon), (
                                gene_name, codon_index, orig_codon,
                                recoded_codon)
            except AssertionError as e:
                if (gene_name in MANUAL_SECIS_FIXES or
                        gene_name in MANUAL_AGR_FIXES):
                    pass  # ok
                elif (orig_codon in START_CODON_SET and
                        recoded_codon in START_CODON_SET):
                    pass  # ok
                else:
                    raise e

            orig_pos_pythonic = get_global_position(
                    feature, maybe_prfB_hack_codon_index,
                    validation_global_seq=orig_seq,
                    validation_sub_seq=orig_codon)

            recoded_pos_pythonic = get_global_position(
                    recoded_f, codon_index, validation_global_seq=recoded_seq,
                    validation_sub_seq=recoded_codon)

            observed_forbidden_codon_positions.add(int(recoded_pos_pythonic))

            data_obj = {
                'segment': map_gene_to_segment[gene_name],
                'gene': gene_name,
                'codon_index': codon_index,
                'strand': feature.strand,
                'position_mds42': orig_pos_pythonic + 1,
                'position_rc3': recoded_pos_pythonic + 1,
                'codon_mds42': orig_codon,
                'codon_rc3': recoded_codon,
            }

            data_obj_list.append(data_obj)

            codon_index += 3

    # Make sure we got all of them.
    missing_forbidden_codon_positions = (
            forbidden_codon_positions - observed_forbidden_codon_positions)
    assert not missing_forbidden_codon_positions

    # DEBUG
    # print ('missing', len(missing_forbidden_codon_positions),
    #         sorted(missing_forbidden_codon_positions))

    df = pd.DataFrame(data_obj_list)
    df.to_csv(output_report, index=False)


def count_nucleotide_changes_in_coding_regions(
        orig_seq_record, recoded_seq_record,
        annotations_fixed_recoded_seq_record):
    """Count number of nucleotide changes that fall in codons.

    NOTE: Partially copied from generate_list_of_recoded_codons()
    """
    orig_seq = str(orig_seq_record.seq)
    recoded_seq = str(recoded_seq_record.seq)

    codon_usage_memex_obj = get_ecoli_codon_usage_memex(randomized=True)

    START_CODON_SET = set(codon_usage_memex.CodonUsageMemex.start_codons)

    data_obj_list = []

    orig_cds_features = [f for f in orig_seq_record.features if f.type == 'CDS']

    # Manually include prfB gene. No CDS in original MDS42 record because of
    # frameshift.
    orig_cds_features.append(
            get_feature_by_gene_name('prfB', 'gene', orig_seq_record))
    # orig_cds_features = [get_feature_by_gene_name('prfB', 'gene', orig_seq_record)]
    orig_cds_features = sorted(
            orig_cds_features, key=lambda f: f.location.start)

    total_coding_features = len(orig_cds_features)

    # Convenience map to make gene lookups below constant time.
    map_gene_name_to_feature_in_recoded = {}
    for f in recoded_seq_record.features:
        if not f.type == 'CDS':
            continue
        gene_name = get_feature_gene_name(f)
        assert not gene_name in map_gene_name_to_feature_in_recoded, gene_name
        map_gene_name_to_feature_in_recoded[gene_name] = f

    # Annotated forbidden codon positions.
    forbidden_codon_positions = set()
    for f in annotations_fixed_recoded_seq_record.features:
        if f.type == 'recoded_codon':
            if f.strand == 1:
                forbidden_codon_positions.add(int(f.location.start))
            else:
                forbidden_codon_positions.add(int(f.location.start) + 2)
    assert len(forbidden_codon_positions) == 62214

    # Track codon positions we've seen and compare to forbiddenCodon
    # annotations.
    observed_forbidden_codon_positions = set()

    total_nt_count = 0
    total_nt_delta = 0
    forbidden_nt_count = 0
    forbidden_nt_delta = 0

    for idx in range(0, total_coding_features):
        feature = orig_cds_features[idx]

        # DEBUG: Verbose output.
        # if idx % 100 == 0:
        #     print 'Analyzing %d-strand feature %s,  %d of %d' % (
        #             feature.strand, get_feature_gene_name(feature),
        #             idx + 1, total_coding_features)

        gene_name = get_feature_gene_name(feature)
        assert gene_name is not None

        recoded_f = map_gene_name_to_feature_in_recoded.get(gene_name)
        if recoded_f is None:
            print 'Gene %s (%d) not found in recoded genome. Skipping.' % (
                    gene_name, feature.location.start)
            continue

        if gene_name == 'prfB':
            assert len(feature) == len(recoded_f) + 1
        else:
            assert len(feature) == len(recoded_f)

        orig_feature_seq = feature.extract(orig_seq)
        recoded_feature_seq = recoded_f.extract(recoded_seq)

        codon_index = 0
        while codon_index < len(feature):

            # Manually frameshift.
            maybe_prfB_hack_codon_index = codon_index
            if gene_name == 'prfB' and codon_index > 72:
                maybe_prfB_hack_codon_index += 1

            orig_codon = orig_feature_seq[
                    maybe_prfB_hack_codon_index:maybe_prfB_hack_codon_index + 3]

            recoded_codon = recoded_feature_seq[codon_index:codon_index + 3]

            total_nt_count += 3
            total_nt_delta += _get_nt_delta(orig_codon, recoded_codon)

            # If forbidden, verify synonymous and update.
            if orig_codon in FORBIDDEN_CODONS:
                try:
                    assert codon_usage_memex_obj.are_synonymous(
                            orig_codon, recoded_codon), (
                                    gene_name, codon_index, orig_codon,
                                    recoded_codon)
                except AssertionError as e:
                    if (gene_name in MANUAL_SECIS_FIXES or
                            gene_name in MANUAL_AGR_FIXES):
                        pass  # ok
                    elif (orig_codon in START_CODON_SET and
                            recoded_codon in START_CODON_SET):
                        pass  # ok
                    else:
                        raise e

                forbidden_nt_count += 3
                forbidden_nt_delta += _get_nt_delta(orig_codon, recoded_codon)

            recoded_pos_pythonic = get_global_position(
                    recoded_f, codon_index, validation_global_seq=recoded_seq,
                    validation_sub_seq=recoded_codon)

            observed_forbidden_codon_positions.add(int(recoded_pos_pythonic))

            codon_index += 3

    # Make sure we got all of them.
    missing_forbidden_codon_positions = (
            forbidden_codon_positions - observed_forbidden_codon_positions)
    assert not missing_forbidden_codon_positions

    # DEBUG
    # print ('missing', len(missing_forbidden_codon_positions),
    #         sorted(missing_forbidden_codon_positions))

    return {
        'total_nt_count': total_nt_count,
        'total_nt_delta': total_nt_delta,
        'forbidden_nt_count': forbidden_nt_count,
        'forbidden_nt_delta': forbidden_nt_delta,
    }


def _get_nt_delta(a, b):
    """Returns number of differences between two sequences.
    """
    assert len(a) == len(b)
    delta = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            delta += 1
    return delta


def make_genbank_with_forbidden_annotated(mds42_seq_record):
    """Generates a Genbank with forbidden codons annotated.
    """
    MDS42_FORBIDDEN_ANNOTATED_GENBANK = os.path.join(
            GENOMES_DIR, 'mds42', 'mds42_forbidden_annotated.gbk')

    mds42_seq_str = str(mds42_seq_record.seq)

    new_forbidden_codon_features = []

    cds_features = [f for f in mds42_seq_record.features if f.type == 'CDS']

    for mds42_f in cds_features:
        mds42_f_seq = mds42_f.extract(mds42_seq_str)
        for codon_start in range(3, len(mds42_f_seq), 3):
            mds42_codon = mds42_f_seq[codon_start:codon_start + 3]
            if not mds42_codon in FORBIDDEN_CODONS:
                continue
            mds42_global_pos = get_global_position(mds42_f, codon_start)
            if mds42_f.strand == 1:
                codon_interval = (mds42_global_pos, mds42_global_pos + 3)
            else:
                codon_interval = (mds42_global_pos - 2, mds42_global_pos + 1)
            new_forbidden_codon_features.append(SeqFeature(
                    location=FeatureLocation(*codon_interval),
                    strand=mds42_f.strand,
                    type='forbidden_codon'
            ))

    mds42_seq_record.features += new_forbidden_codon_features
    mds42_seq_record.features = sorted(
            mds42_seq_record.features, key=lambda f: int(f.location.start))

    with open(MDS42_FORBIDDEN_ANNOTATED_GENBANK, 'w') as output_fh:
        SeqIO.write(mds42_seq_record, output_fh, 'genbank')


if __name__ == '__main__':
    MDS42_GENBANK = os.path.join(DATA_DIR, 'mds42_full.gbk')

    print 'Reading MDS42 genbank...'
    mds42_seq_record = SeqIO.read(MDS42_GENBANK, 'gb')
    print '...Done.'

    make_genbank_with_forbidden_annotated(mds42_seq_record)

