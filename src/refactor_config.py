"""
Central place for storing parameters for a refactor.

NOTE: Previously we were using the refactor_config module for this purpose
but it has grown unwieldy at this point. So for now we have a bit of a messy
hybrid.

NOTE(8/26/13): We are now transitioning toward using .yaml config files, as
opposed to this mess. In the interim we'll have an hybrid solution as we
move more parameters to the config file.
"""

import os
import pickle
import pprint
import yaml

from Bio import Restriction
from Bio.Seq import reverse_complement
from Bio.SeqFeature import FeatureLocation

from biopython_util import get_genome_record
from biopython_util import get_feature_by_id
from biopython_util import swap_feature_codon_at_position
from codon_usage_memex import build_codon_usage_dict
from codon_usage_memex import CodonUsageMemex
from manual_modifications import handle_bulk_manual_swaps
from manual_modifications import handle_AGR_problem_ids

from paths import CACHE_DIR
from paths import CONFIG_DIR
from paths import DATA_DIR
from paths import EC_ALIGNMENT_DATA_DIR
from paths import GENOMES_DIR
from paths import OUTPUT_DIR
from paths import PWD
from paths import TMP_DATA_DIR


# Make sure NUPACKHOME environment variable is set.
# This is required by the Salis rbs_calc package.
os.environ["NUPACKHOME"] = "/opt/nupack"

###############################################################################
# The yaml config we are using.
# TODO: Transition to make this module into a config model wrapper which
#     is instantiated from main.
###############################################################################

CURRENT_YAML_CONFIG = os.path.join(CONFIG_DIR, 'seven_codon_config.yaml')
with open(CURRENT_YAML_CONFIG) as yaml_config_fh:
    YAML_CONFIG_DICT = yaml.load(yaml_config_fh)


###############################################################################
# Misc (to be organized)
###############################################################################

DEBUG = False

DEBUG_SINGLE_ITERATION = False

PRETTY_PRINTER = pprint.PrettyPrinter(indent=4)

USE_CACHE = True

GENOME_SOURCE = os.path.join(DATA_DIR, 'mds42_full.gbk')

NUM_CORES = YAML_CONFIG_DICT.get('num_cores', 1)

SELENOCYSTEINE_CODON = 'TGA'

# Design constraints based on Gen9 Rules.
# See http://gen9bio.com/faq/
HOMOPOLYMER_RUN_LIMIT_A = 8
HOMOPOLYMER_RUN_LIMIT_C = 8
HOMOPOLYMER_RUN_LIMIT_G = 5
HOMOPOLYMER_RUN_LIMIT_T = 8

CODONS_TO_REMOVE = YAML_CONFIG_DICT['forbidden_codons']
CODON_USAGE_REPORT = os.path.join(CONFIG_DIR, YAML_CONFIG_DICT['codon_usage'])

ORIGINAL_CODON_USAGE_MEMEX = CodonUsageMemex(build_codon_usage_dict())

REFACTORED_CODON_USAGE_MEMEX = (
        CodonUsageMemex.build_from_removed_codons_usage_report(
                CODON_USAGE_REPORT))


###############################################################################
# Manual feature handling.
#
# Functions that capture manual changes that aren't handled by any other part
# of the pipeline. These should be executed before the rest of the refactoring
# pipeline. There might be a more generic way to do this, but just doing it
# manually for now until we figure that out.
###############################################################################

def handle_prfB(genome_record):
    """Modifies the genome record to make prfB a CDS and remove
    the frameshift.
    """
    prfB_locus_tag = 'ECMDS42_2390'
    prfB_feature = get_feature_by_id(genome_record, prfB_locus_tag)
    prfB_feature_seq = str(prfB_feature.extract(genome_record.seq))
    original_length = len(prfB_feature)

    # Feature is on the negative strand.
    assert prfB_feature.strand == -1

    # First delete the frameshift, CTTT -> CTT, by removing the T at the 75th
    # position of the feature.
    FRAMESHIFT_START = 72
    assert 'CTTTGAC' == str(prfB_feature_seq[FRAMESHIFT_START:FRAMESHIFT_START + 7])
    prfB_feature_seq = prfB_feature_seq[:74] + prfB_feature_seq[75:]
    assert 'CTTGAC' == str(prfB_feature_seq[FRAMESHIFT_START:FRAMESHIFT_START + 6])

    # Next change the RBS just before it AGGGGG -> CGTGGG (as per Chris Gregg
    # data from 8/25/13).
    RBS_START = 63
    assert 'AGGGGG' == str(prfB_feature_seq[63:69])
    prfB_feature_seq = prfB_feature_seq[:63] + 'CGT' + prfB_feature_seq[66:]
    assert 'CGTGGG' == str(prfB_feature_seq[63:69])

    # Now replace the underlying sequence of the genome.
    updated_seq = (
            genome_record.seq[:prfB_feature.location.start] +
            reverse_complement(prfB_feature_seq) +
            genome_record.seq[prfB_feature.location.end:]
    )
    assert len(genome_record.seq) - 1 == len(updated_seq)
    genome_record.seq = updated_seq

    # Change the prfB feature to be a CDS.
    prfB_feature.type = 'CDS'

    # Bump the start position one unit right.
    # Remember, it's on the negative strand.
    prfB_feature.location = FeatureLocation(
            prfB_feature.location.start,
            prfB_feature.location.end - 1,
            strand=prfB_feature.strand
    )

    # Update the positions of downstream features.
    updated_features = []
    for feature in genome_record.features:
        if feature.location.start > prfB_feature.location.start:
            feature = feature._shift(-1)
        updated_features.append(feature)
    genome_record.features = updated_features

    # Make sure changes went through.
    mod_prfB_feature = get_feature_by_id(genome_record, prfB_locus_tag)
    mod_prfB_feature_seq = str(mod_prfB_feature.extract(genome_record.seq))
    assert original_length - 1 == len(mod_prfB_feature)
    assert prfB_feature_seq[:6] == mod_prfB_feature_seq[:6], (
            "Before: %s, After: %s" %
            (prfB_feature_seq[:6], mod_prfB_feature_seq[:6]))
    assert prfB_feature_seq[-10:] == mod_prfB_feature_seq[-10:], (
            "Before: %s, After: %s" %
            (prfB_feature_seq[-10:], mod_prfB_feature_seq[10:]))

    return genome_record


def handle_fdhF(genome_record):
    """Perform amino acid changes to preserve the SECIS site sufficiently.

    Gly -> Gly (C424T)
    Val -> Ile (G454A)

    Returns the updated genome_record.
    """
    feature_id = 'ECMDS42_3518'

    CODON_POS = 423
    PREVIOUS_CODON = 'GGC'
    NEW_CODON = 'GGT'
    genome_record = swap_feature_codon_at_position(
            genome_record, feature_id, CODON_POS, PREVIOUS_CODON, NEW_CODON)

    CODON_POS = 453
    PREVIOUS_CODON = 'GTC'
    NEW_CODON = 'ATC'
    return swap_feature_codon_at_position(
            genome_record, feature_id, CODON_POS, PREVIOUS_CODON, NEW_CODON)


def handle_fdnG(genome_record):
    """Perform amino acid changes to preserve the SECIS site sufficiently.

    Ser -> Ala

    Returns the updated genome_record.
    """
    feature_id = 'ECMDS42_1186'

    CODON_POS = 606
    PREVIOUS_CODON = 'AGT'
    NEW_CODON = 'GCT'
    return swap_feature_codon_at_position(
            genome_record, feature_id, CODON_POS, PREVIOUS_CODON, NEW_CODON)


def handle_fdoG(genome_record):
    """Perform amino acid changes to preserve the SECIS site sufficiently.

    Ser -> Ala

    Returns the updated genome_record.
    """
    feature_id = 'ECMDS42_3333'

    CODON_POS = 606
    PREVIOUS_CODON = 'AGT'
    NEW_CODON = 'GCT'
    return swap_feature_codon_at_position(
            genome_record, feature_id, CODON_POS, PREVIOUS_CODON, NEW_CODON)


MANUAL_CHANGE_FUNCTIONS = [
        handle_prfB,
        handle_fdhF,
        handle_fdnG,
        handle_fdoG
]

###############################################################################
# Original genome record
#
# The original genome record is the one read from source, with the manual
# update functions run over it.
###############################################################################

def create_original_genome_record():
    """Creates the original genome record.
    """
    genome_record = get_genome_record(GENOME_SOURCE)
    for fn in MANUAL_CHANGE_FUNCTIONS:
        genome_record = fn(genome_record)

    bulk_manual_swaps_file = YAML_CONFIG_DICT.get('bulk_manual_fixes', None)
    if bulk_manual_swaps_file:
        MG1655_SOURCE = os.path.join(GENOMES_DIR, 'mg1655', 'mg1655.genbank')
        MG1655_GENOME_RECORD = get_genome_record(MG1655_SOURCE,
                use_old_id_strategy=True)
        bulk_manual_swaps_full_path = os.path.join(
                DATA_DIR, YAML_CONFIG_DICT['bulk_manual_fixes'])
        handle_bulk_manual_swaps(genome_record, bulk_manual_swaps_full_path,
                MG1655_GENOME_RECORD)
        # HACK: Do the tricky ones manually. There's only five of them.
        handle_AGR_problem_ids(genome_record)
    return genome_record


ORIGINAL_GENOME_RECORD = create_original_genome_record()

# HACK
TRANSLATION_CONSERVED_EXCEPTIONS = set([
        'ECMDS42_1886', # folC, fixed in AGR experiments
        'ECMDS42_3717', # not sure why this is breaking but it's outside our synthesis bounds
])


###############################################################################
# Restriction enzymes to remove.
###############################################################################

RESTRICTION_ENZYME_SITES_TO_REMOVE = [
        Restriction.BsaI,
        Restriction.BsmBI,
        Restriction.AarI
]

###############################################################################
# Manual edits
###############################################################################

AGN_SEPARATION_DATA_FILE = os.path.join(CONFIG_DIR, 'AGN_debug_MJL.csv')
