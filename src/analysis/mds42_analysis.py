"""
Script to analyze overlaps in MDS42.
"""

import csv
import os

import yaml

from analysis import determine_codon_usage_table
from analysis import determine_dicodon_usage_dict
from analysis import determine_feature_frame
from analysis import determine_hexamer_usage_dict
from biopython_util import get_feature_gene
from biopython_util import get_genome_record
from codon_usage_memex import CodonUsageMemex
from conflicting_pair_finder import find_all_overlaps
from paths import CONFIG_DIR
from paths import DATA_DIR
from tool_wrappers.rbs_calc_util import calc_rbs_score_for_feature


def analyze_conflicting_pairs():
    genome_record = get_genome_record('../data/mds42_full.gbk')

    # Write the results to a file.
    OUTPUT_FILE = '../data/mds42_conflicting_pair_analysis.csv'
    FIELD_NAMES = [
            'conflict_type',
            'overlap_type',
            'overlap_size',
            'same_start',
            'same_end',
            'polarity',
            'upstream_type',
            'upstream_locus',
            'upstream_gene',
            'upstream_start',
            'upstream_end',
            'upstream_frame',
            'downstream_type',
            'downstream_locus',
            'downstream_gene',
            'downstream_start',
            'downstream_end',
            'downstream_frame',
    ]

    with open(OUTPUT_FILE, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, FIELD_NAMES)
        writer.writeheader()

        # Figure out all the overlaps.
        overlap_list = find_all_overlaps(
                genome_record,
                CODONS_TO_REMOVE,
                require_feature_type=[],
                ignore_feature_type=['source', 'gene', 'misc_feature'],
                essential_feature_ids=[],
                full_analysis=True,
                cache=False,
                include_close_features=True,
        )

        # Write the results.
        for overlap_obj in overlap_list:
            upstream_feature = overlap_obj['upstream_feature']
            downstream_feature = overlap_obj['downstream_feature']
            conflict_type = overlap_obj['conflict_type']
            overlap_type = overlap_obj.get('overlap_type', 'NA')
            overlap_size = overlap_obj.get('overlap_size', 'NA')
            polarity = overlap_obj.get('polarity', 'NA')

            upstream_gene = get_feature_gene(upstream_feature)
            if not upstream_gene:
                upstream_gene = 'UNKNOWN'

            downstream_gene = get_feature_gene(downstream_feature)
            if not downstream_gene:
                downstream_gene = 'UNKNOWN'

            same_start = (
                    upstream_feature.location.start ==
                            downstream_feature.location.start)

            same_end = (
                    upstream_feature.location.end ==
                            downstream_feature.location.end)

            row = {
                    'conflict_type': conflict_type,
                    'overlap_type': overlap_type,
                    'overlap_size': overlap_size,
                    'polarity': polarity,
                    'same_start': same_start,
                    'same_end': same_end,
                    'upstream_type': upstream_feature.type,
                    'upstream_locus':
                            upstream_feature.qualifiers.get(
                                    'locus_tag', ['UNKNOWN'])[0],
                    'upstream_gene': upstream_gene,
                    'upstream_start': upstream_feature.location.start,
                    'upstream_frame':
                            determine_feature_frame(upstream_feature),
                    'upstream_end': upstream_feature.location.end,
                    'downstream_type': downstream_feature.type,
                    'downstream_locus':
                            downstream_feature.qualifiers.get(
                                    'locus_tag', ['UNKNOWN'])[0],
                    'downstream_gene': downstream_gene,
                    'downstream_start': downstream_feature.location.start,
                    'downstream_end': downstream_feature.location.end,
                    'downstream_frame':
                            determine_feature_frame(downstream_feature),
            }
            writer.writerow(row)


def analyze_codon_usage(refactor_config_yaml, original_record_path,
        refactored_record_path):
    """Analyze MDS42 codon usage before and after and refactoring.

    Args:
        refactor_config_yaml: File containing the configuration for this
            particular refactor (e.g. codons to remove, etc.).
    """
    with open(refactor_config_yaml) as yaml_fh:
        YAML_CONFIG_DICT = yaml.load(yaml_fh)

    CODONS_TO_REMOVE = YAML_CONFIG_DICT['forbidden_codons']

    ignore_feature_ids = [
            'CDS_fdnG_ECMDS42_1186_1254228_1257276',
            'CDS_fdhF_ECMDS42_3518_3726321_3728469',
            'CDS_fdoG_ECMDS42_3333_3511985_3515036'
    ]

    original_genome_record = get_genome_record(original_record_path,
            use_old_id_strategy=True, only_get_features=['CDS'])

    refactored_genome_record = get_genome_record(refactored_record_path,
            use_old_id_strategy=True, only_get_features=['CDS'])

    OUTPUT_FILE = refactored_record_path + '.codon_usage_analysis.csv'
    FIELD_NAMES = [
            'codon',
            'attempted_remove',
            'amino_acid',
            'original_usage',
            'refactored_usage',
            'target_usage',
            'delta_usage',
            'original_count',
            'original_amino_acid_count',
            'refactored_count',
            'refactored_amino_acid_count',
    ]

    CODON_USAGE_REPORT = os.path.join(CONFIG_DIR,
            YAML_CONFIG_DICT['codon_usage'])
    TARGET_CODON_USAGE_MEMEX = (
            CodonUsageMemex.build_from_removed_codons_usage_report(
                    CODON_USAGE_REPORT))

    print 'Calculating original source usage...'
    original_codon_usage = determine_codon_usage_table(
            original_genome_record)

    print 'Calculating refactored genome codon usage...'
    refactored_codon_usage = determine_codon_usage_table(
            refactored_genome_record,
            ignore_feature_ids=ignore_feature_ids)

    print 'Writing result ...'
    with open(OUTPUT_FILE, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, FIELD_NAMES)
        writer.writeheader()

        for codon, original_data in original_codon_usage.iteritems():
            refactored_data = refactored_codon_usage[codon]
            row = {
                'codon': codon,
                'attempted_remove': codon in CODONS_TO_REMOVE,
                'amino_acid': original_data['amino_acid'],
                'original_usage': "{0:.2f}".format(original_data['usage']),
                'refactored_usage': "{0:.2f}".format(refactored_data['usage']),
                'target_usage':
                        TARGET_CODON_USAGE_MEMEX.get_codon_usage(codon),
                'delta_usage': "{0:.2f}".format(
                        refactored_data['usage'] - original_data['usage']),
                'original_count': original_data['count'],
                'original_amino_acid_count': original_data['amino_acid_count'],
                'refactored_count': refactored_data['count'],
                'refactored_amino_acid_count': refactored_data['amino_acid_count'],
            }
            writer.writerow(row)


def analyze_hexamer_usage():
    """Analyze MDS42 hexamer usage before and after and refactoring.
    """
    print 'Analyzing hexamer usage ...'

    MDS42_GENOME_SOURCE_ORIGINAL = '../data/mds42_full.gbk'
    original_genome_record = get_genome_record(MDS42_GENOME_SOURCE_ORIGINAL)

    MDS42_GENOME_SOURCE_REFACTORED = (
            '../data/2013_01_28_20_40_42_mds42_refactored.gbk')
    refactored_genome_record = get_genome_record(MDS42_GENOME_SOURCE_REFACTORED)

    OUTPUT_FILE = '../data/mds42_hexamer_usage_analysis.csv'
    FIELD_NAMES = [
            'hexamer',
            'original_usage',
            'refactored_usage',
            'delta_usage',
            'original_count',
            'refactored_count',
    ]

    with open(OUTPUT_FILE, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, FIELD_NAMES)
        writer.writeheader()

        original_hexamer_usage = determine_hexamer_usage_dict(
                original_genome_record)
        refactored_hexamer_usage = determine_hexamer_usage_dict(
                refactored_genome_record)

        NO_USAGE_DICT = {
                'usage': 0.0,
                'count': 0,
        }

        all_hexamers = (set(original_hexamer_usage.keys()) |
                set(refactored_hexamer_usage.keys()))

        for hexamer in all_hexamers:
            original_data = original_hexamer_usage.get(
                    hexamer, NO_USAGE_DICT)
            refactored_data = refactored_hexamer_usage.get(
                    hexamer, NO_USAGE_DICT)
            row = {
                'hexamer': hexamer,
                'original_usage': original_data['usage'],
                'refactored_usage': refactored_data['usage'],
                'delta_usage': (
                        refactored_data['usage'] - original_data['usage']),
                'original_count': original_data['count'],
                'refactored_count': refactored_data['count'],
            }
            writer.writerow(row)

    print 'Done.'


def analyze_dicodon_usage():
    """Analyze MDS42 dicodon usage before and after and refactoring.

    NOTE: Probably should reuse code from analyze_hexamer_usage().
    """
    print 'Analyzing di-codon usage ...'

    MDS42_GENOME_SOURCE_ORIGINAL = '../data/mds42_full.gbk'
    original_genome_record = get_genome_record(MDS42_GENOME_SOURCE_ORIGINAL)

    MDS42_GENOME_SOURCE_REFACTORED = (
            '../data/2013_01_28_20_40_42_mds42_refactored.gbk')
    refactored_genome_record = get_genome_record(MDS42_GENOME_SOURCE_REFACTORED)

    OUTPUT_FILE = '../data/mds42_dicodon_usage_analysis.csv'
    FIELD_NAMES = [
            'dicodon',
            'original_usage',
            'refactored_usage',
            'delta_usage',
            'original_count',
            'refactored_count',
    ]

    with open(OUTPUT_FILE, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, FIELD_NAMES)
        writer.writeheader()

        original_hexamer_usage = determine_dicodon_usage_dict(
                original_genome_record)
        refactored_hexamer_usage = determine_dicodon_usage_dict(
                refactored_genome_record)

        NO_USAGE_DICT = {
                'usage': 0.0,
                'count': 0,
        }

        all_hexamers = (set(original_hexamer_usage.keys()) |
                set(refactored_hexamer_usage.keys()))

        for hexamer in all_hexamers:
            original_data = original_hexamer_usage.get(
                    hexamer, NO_USAGE_DICT)
            refactored_data = refactored_hexamer_usage.get(
                    hexamer, NO_USAGE_DICT)
            row = {
                'dicodon': hexamer,
                'original_usage': original_data['usage'],
                'refactored_usage': refactored_data['usage'],
                'delta_usage': (
                        refactored_data['usage'] - original_data['usage']),
                'original_count': original_data['count'],
                'refactored_count': refactored_data['count'],
            }
            writer.writerow(row)

    print 'Done.'


def analyze_codons_removed():
    """Analyze what codons were actually removed.

    This is in response to finding that the codon usage following refactoring
    doesn't seem to reflect the removed codons very well.
    """
    print 'Analyzing codons removed...'

    CODONS_TO_REMOVE = [
        'ACC', # T
        'AGA', # R
        'AGG', # R
        'AGT', # S'
        'AGC', # S'
        'ATA', # I
        'CCC', # P
        'CGG', # R
        'GCC', # A
        'GTC', # V
        'TAG', # B'
        'TGA', # U
        'TTA', # L'
        'TTG', # L'
        'TCC', # S
    ]

    MDS42_GENOME_SOURCE_REFACTORED = (
            '../data/2013_01_28_20_40_42_mds42_refactored.gbk')
    refactored_genome_record = get_genome_record(MDS42_GENOME_SOURCE_REFACTORED)

    # Codons only removed from CDS features.
    coding_features = filter(lambda feature: feature.type == 'CDS',
            refactored_genome_record.features)
    print 'Total CDS features', len(coding_features)

    START_CODONS = ['ATG', 'GTG', 'GTC', 'GTT']
    for feature in coding_features:
        seq = feature.extract(refactored_genome_record).seq
        start_codon = str(seq[0:3])
        if start_codon not in START_CODONS:
            print 'strand', feature.strand
            print start_codon
            print feature
            print seq
            print '\n\n'
            print str(refactored_genome_record[
                    feature.location.start:feature.location.end].seq)
            break


def analyze_rbs_strengths():
    """Analyze RBS strengths in the original MDS42 genome.
    """
    genome_record = get_genome_record('../data/mds42_full.gbk')

    forward_seq = genome_record.seq

    coding_features = filter(
            lambda f: f.type == 'CDS',
            genome_record.features)

    # Get RBS scores for entire forward and reverse sequence.

    RBS_CALC_BUFFER = 30

    OUTPUT_FILE = '../data/mds42_rbs_strength_analysis_size_30_buffer.csv'
    FIELD_NAMES = [
            'locus_tag',
            'gene',
            'rbs',
            'dG',
            'kinetic_score',
            'expression',
            'seq'
    ]

    with open(OUTPUT_FILE, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, FIELD_NAMES)
        writer.writeheader()

        num_features = len(coding_features)
        for feature_index in range(num_features):
            print 'Calculating feature %d of num_features %d' % (
                    feature_index + 1, num_features)

            feature = coding_features[feature_index]
            row = calc_rbs_score_for_feature(
                    feature, forward_seq, RBS_CALC_BUFFER)
            writer.writerow(row)

    print 'Done'


if __name__ == '__main__':
    from paths import OUTPUT_DIR

    CONFIG = os.path.join(CONFIG_DIR, 'seven_codon_config.yaml')
    MDS42_GENOME_SOURCE_ORIGINAL = os.path.join(DATA_DIR, 'mds42_full.gbk')

    SEG5_2MB_DIR = os.path.join(DATA_DIR, 'completed_segments', 'seg5-seg48')
    RECODED_PATH = os.path.join(SEG5_2MB_DIR,
            '2013_08_30.seg5-seg48.motifs_removed.AGN_fixed.genbank')
    analyze_codon_usage(CONFIG, MDS42_GENOME_SOURCE_ORIGINAL, RECODED_PATH)
