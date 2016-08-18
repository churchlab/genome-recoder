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
Methods for performing manual or hard-coded modifications.
"""

import csv
import re

from Bio.Seq import reverse_complement

from biopython_util import get_feature_gene
from biopython_util import insert_sequence_and_update_features
from biopython_util import swap_feature_codon_at_position
from biopython_util import update_feature_seq_at_position


# These are ids of records that I've identified have problems and need
# fixing.
KNOWN_PROBLEM_IDS = set([
    'holB_AGA_4_ATAATGCGC',
    'folC_AGR_16-19_CGGCGA_opt',
    'ssb_AGA_10_CGA_wobble',
    'dnaT_AGA_529_CGC_opt',
    'prfB_FS+CGN+Wobble_opt'
])


def handle_bulk_manual_swaps(genome_record, input_file, mg1655_genome_record):
    """Method that allows handling bulk swaps.

    NOTE: This method is pretty much hard-coded to work with derivatives of
    MG1655.
    """
    # After considering all the various options for finding the exact position
    # to change, I've decided to go with parsing the 'AGR ID', which has the
    # name of the gene as well as the codon position. Let's see how it goes.

    # Eventually we'll need a way to get from gene name to feature in the
    # MDS42 genome.  One problem that arises when matching gene names is we
    # have gene name synonyms.  Thus we can use the MG1655 record which
    # has a lot more synonyms recorded.

    # Create a bi-map linking various gene synonyms of the MG1655 genbank.
    # This allows us to handle more flexible cases.
    # TODO: In general, this synonym-finding functionality could be useful
    # elsewhere. Maybe use Ecocyc or regulondb data for this purpose.
    gene_to_synonym_bimap = {}
    mg1655_cds_features = [feature for feature in mg1655_genome_record.features
            if feature.type in set(['CDS', 'gene'])]
    for feature in mg1655_cds_features:
        if not feature.type in ['CDS', 'gene']:
            continue

        maybe_gene = get_feature_gene(feature)
        if not maybe_gene:
            continue

        if 'gene_synonym' in feature.qualifiers:
            # Build a list containing all synonyms which will serve as the
            # value of the bimap.
            synonym_list = [maybe_gene]

            # Check each to see if it can be split.
            for synonym_phrase in feature.qualifiers['gene_synonym']:
                split_phrase = synonym_phrase.split(';')
                for synonym in split_phrase:
                    clean_synonym = synonym.strip()
                    if len(clean_synonym) > 0:
                        synonym_list.append(clean_synonym)

            for synonym in synonym_list:
                gene_to_synonym_bimap[synonym] = synonym_list

    # Create a map from gene name to CDS feature for that gene in genome_record.
    gene_to_feature_map = {}
    cds_features = [feature for feature in genome_record.features
            if feature.type == 'CDS']
    for feature in cds_features:
        maybe_gene = get_feature_gene(feature)
        if not maybe_gene:
            continue

        # Always add the feature for the gene. It's possible one of the
        # synonyms was added earlier, so we override it.
        gene_to_feature_map[maybe_gene] = feature

        # Get all synonyms and build up the map.
        if 'gene_synonym' in feature.qualifiers:
            synonym_set = set(
                    [maybe_gene] +
                    feature.qualifiers['gene_synonym'] +
                    gene_to_synonym_bimap.get(maybe_gene, []))
            for synonym in synonym_set:
                # Don't override if it's already there. We want actual genes
                # to get precedence over less reliable synonyms.
                if not synonym in gene_to_feature_map:
                    gene_to_feature_map[synonym] = feature

    # And use the bimap to pickup any missing synonym connections.
    for gene, synonym_list in gene_to_synonym_bimap.iteritems():
        if gene in gene_to_feature_map:
            continue
        for synonym in synonym_list:
            if synonym in gene_to_feature_map:
                feature = gene_to_feature_map[synonym]
                gene_to_feature_map[gene] = feature
                break

    # Now iterate through the manual fixes and make the changes.
    with open(input_file) as input_fh:
        reader = csv.DictReader(input_fh, delimiter='\t')
        for manual_fix in reader:
            clean_id = manual_fix['AGR ID'].strip()
            if clean_id in KNOWN_PROBLEM_IDS:
                continue

            parsed_id = re.match(
                    r'(?P<gene>[a-zA-Z]+)_[a-zA-Z]+_(?P<mutation_start>[0-9]+)?.*',
                    clean_id)

            # prfB is the only know weird case, so check that anything weird
            # is indeed prfB and skip it for now.
            if not parsed_id:
                parsed_id = re.match(r'(?P<gene>[a-zA-Z]+)_.*',
                        manual_fix['AGR ID'])
            gene = parsed_id.group('gene')

            # Make it 0-indexed per the method that does the swap's API.
            mutation_start = int(parsed_id.group('mutation_start')) - 1

            if not gene in gene_to_feature_map:
                print "%s not in map." % gene
                assert False

            feature = gene_to_feature_map[gene]

            # Parse the old and new sequence data.
            previous_seq = manual_fix['wt Genotype'].upper()
            new_seq = manual_fix['Destination Genotype'].upper()

            if len(new_seq) == 0:
                continue

            assert len(previous_seq) == 3
            assert len(new_seq) == 3

            # Ugh, looks like the data gives changes always in the forward
            # strand.  We'll have to reverse it ourselves in the negative
            # strand case.
            if feature.strand == -1:
                mutation_start = len(feature) - mutation_start - len(previous_seq)
                previous_seq = reverse_complement(previous_seq)
                new_seq = reverse_complement(new_seq)

            swap_feature_codon_at_position(genome_record, feature.id,
                    mutation_start, previous_seq, new_seq)


def handle_AGR_problem_ids(genome_record):
    """HACK: Manually handle problematic AGR-related swaps for now.
    """
    ### 'prfB_FS+CGN+Wobble_opt' - handled in manual fixes

    ### 'holB_AGA_4_ATAATGCGC'
    # Shift holB
    insert_sequence_and_update_features(genome_record, 'ATA',
            987133, no_feature_shift_if_inside=True)

    # NOTE: Those following holB have position +3 due to insertion.

    ### 'folC_AGR_16-19_CGGCGA_opt'
    update_feature_seq_at_position(genome_record, 2002038 + 3,
            'ACTTACTTGCCACCGCTTCTCCTCGCG',
            'ACTTACTTGCCATCGCTTCGCCGCGCG')

    # NOTE: Those following prfB have position 1 less due to the deletion.

    ### 'ssb_AGA_10_CGA_wobble',
    update_feature_seq_at_position(genome_record, 3703337 - 1 + 3,
            'ATGGCCAGCagaGGCGTAAACAAGGTT',
            'ATGGCGAGTCGAGGTGTTAATAAGGTA')

    ### 'dnaT_AGA_529_CGC_opt',
    update_feature_seq_at_position(genome_record, 3936042 - 1 + 3,
            'CAAAACTCTGGAA',
            'TAGAACGCGTGAG') # NOTE change from source to remove TTA




if __name__ == '__main__':
    import os

    from Bio import SeqIO
    import yaml

    from biopython_util import get_genome_record
    from paths import CONFIG_DIR
    from paths import DATA_DIR
    from paths import GENOMES_DIR

    # CURRENT_YAML_CONFIG = os.path.join(CONFIG_DIR, 'seven_codon_config.yaml')
    # with open(CURRENT_YAML_CONFIG) as yaml_config_fh:
    #     YAML_CONFIG_DICT = yaml.load(yaml_config_fh)

    # GENOME_SOURCE = os.path.join(DATA_DIR, 'mds42_full.gbk')
    # genome_record = get_genome_record(GENOME_SOURCE)

    # MG1655_SOURCE = os.path.join(GENOMES_DIR, 'mg1655', 'mg1655.genbank')
    # mg1655_genome_record = get_genome_record(MG1655_SOURCE,
    #         use_old_id_strategy=True)

    # bulk_manual_swaps_file = YAML_CONFIG_DICT.get('bulk_manual_fixes', None)
    # bulk_manual_swaps_full_path = os.path.join(
    #         DATA_DIR, YAML_CONFIG_DICT['bulk_manual_fixes'])

    # handle_bulk_manual_swaps(genome_record, bulk_manual_swaps_full_path,
    #         mg1655_genome_record)

    # handle_AGR_problem_ids(genome_record)
    from refactor_config import ORIGINAL_GENOME_RECORD

    TEST_OUTPUT = os.path.join(DATA_DIR, 'test_manual_changes.gbk')
    with open(TEST_OUTPUT, 'w') as fh:
        SeqIO.write(ORIGINAL_GENOME_RECORD, TEST_OUTPUT, 'genbank')
