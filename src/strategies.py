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
Strategies for refactoring.

The end-to-end refactoring flow will be captured in the methods here.
"""

import copy
import os
import pickle
import re
import stat

from Bio import SeqIO

from biopython_util import get_genome_record
from biopython_util import InsertType
from codon_replacer_util import replace_forbidden_codons
from conflicting_pair_fixer import ConflictingPairFixer
from homology import fix_homology_issues
from refactor_checker import check_all
from refactor_checker import check_forbidden_codons_removed
from refactor_config import AGN_SEPARATION_DATA_FILE
from refactor_config import CODONS_TO_REMOVE
from refactor_config import DEBUG_SINGLE_ITERATION
from refactor_config import NUM_CORES
from refactor_config import ORIGINAL_CODON_USAGE_MEMEX
from refactor_config import OUTPUT_DIR
from refactor_config import REFACTORED_CODON_USAGE_MEMEX
from refactor_context import RefactorContext
from refactor_config import USE_CACHE


def refactor_with_min_overlap_fixes_and_preserve_rbs(
        genome_record,
        tmp_file_prefix,
        debug=False):
    """Refactoring strategy that only fixes overlaps when it's necessary
    for removing forbidden codons and/or preserving rbs strength.

    This is the second big iteration of our strategy following disucssion
    on 2/6/13. Notable differences from the initial strategy:
        * Before this we were just pulling apart all overlaps and copying the
          RBS site. However, intuitively we have concerns that this may
          introduce many unnecessary changes.
        * We are now also taking into account coding features that are close
          enough to each other (even if not overlapping) where re-coding
          could affect one of the feature's RBS regions.

    Super high-level algorithm overview:
        1. Fix overlaps.
        2. Recode each gene to remove forbidden codons.

    Medium-level algorithm overview:
        * Identify all pairs of coding regions that are either overlapping,
          or are close enough (< 20 bp), where recoding one may affect
          translation of the other.
            - For each of these pairs:
                * If there are no forbidden codons in the affected regions:
                    - Nothing to do, mark the pair as resolved.
                * If there are forbidden codons:
                    - Do an exhaustive search over the affected region
                      and try to find a path of synonymous codon
                      substitutions that don't require physical separation.
                      If success, perform the change, and mark any changed
                      codons as "fixed" so that they are not changed in the
                      second half of the overall algorithm where we do the bulk
                      forbidden codon removal.
                    - Otherwise, we need to separate:
                        * If overlap < 4 bp, find minimum amount to copy that
                          resolves any issues, and lock affected RBS regions
                          so that they are not changed.
                        * Otherwise, need to copy overlap + 15 bp upstream
                          of ATG, and to help prevent snap-back:
                            - muddle old start codon in upstream gene
                            - muddle bases in copied region that are not part
                              of RBS
        * Perform synonymous swapping as before, but this time respecting
          locked-in regions from first half of algorithm.
    """
    # Make a copy of the original for validation.
    original_genome_record = copy.deepcopy(genome_record)

    # Context object to be passed around to different methods.
    refactor_context = RefactorContext(genome_record)

    ###########################################################################
    # Fix overlaps
    ###########################################################################

    cfp = ConflictingPairFixer(
            refactor_context,
            cache=USE_CACHE,
            include_close_features=True,
            single_iteration=DEBUG_SINGLE_ITERATION,
            force_separate_AGN=True,
            agn_separation_data_file=AGN_SEPARATION_DATA_FILE
    )
    genome_record = cfp.fix_overlaps()
    refactor_context.set_genome_record(genome_record)

    # Write the output before going on to next step in case there is an error
    # later, so we can at least have partial results.
    _write_output(genome_record, {}, tmp_file_prefix)

    # Validate that overlaps were fixed correctly before moving on.
    check_all(original_genome_record, genome_record)


    ###########################################################################
    # Swap out remaining forbidden codons
    # NOTE: Some were already replaced while fixing overlaps.
    ###########################################################################

    (genome_record, metadata) = replace_forbidden_codons(
            refactor_context,
            num_cores=NUM_CORES,
            tmp_file_prefix=tmp_file_prefix,
            debug=debug)

    # Write the output just in case again.
    _write_output(genome_record, metadata, tmp_file_prefix)

    # Check that forbidden codons were removed.
    check_forbidden_codons_removed(genome_record, CODONS_TO_REMOVE)

    # Validate that we're still good after codon replacement.
    check_all(original_genome_record, genome_record)


    ###########################################################################
    # Resolve homology issues
    ###########################################################################

    genome_record = fix_homology_issues(genome_record)

    # Write the final output, overriding the intermediate write above.
    _write_output(genome_record, metadata, tmp_file_prefix)

    # Validation checks.
    check_forbidden_codons_removed(genome_record, CODONS_TO_REMOVE)
    check_all(original_genome_record, genome_record)

    print 'Done.'


def _write_output(genome_record, metadata, tmp_file_prefix):
    """Write the output to a local file."""

    # Make sure the output directory exists.
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
        # User and group have all permissions | get group id from directory.
        os.chmod(OUTPUT_DIR, stat.S_ISGID | 0775)

    # Write the resulting genome.
    genome_output_file = os.path.join(
            OUTPUT_DIR, tmp_file_prefix + '_mds42_refactored.gbk')
    with open(genome_output_file, 'w') as output_fh:
        SeqIO.write(genome_record, output_fh, 'genbank')

    # Write the metadata.
    metadata_output_file = os.path.join(
            OUTPUT_DIR, tmp_file_prefix + '_mds42_refactored.meta')
    with open(metadata_output_file, 'w') as metadata_output_fh:
        pickle.dump(metadata, metadata_output_fh)


def _get_genome_and_metadata_from_output(tmp_file_prefix):
    """Utility method for getting the genome_record and metadata from
    intermediate output.
    """
    genome_output_file = os.path.join(
            OUTPUT_DIR, tmp_file_prefix + '_mds42_refactored.gbk')
    genome_record = get_genome_record(genome_output_file)

    metadata_output_file = os.path.join(
            OUTPUT_DIR, tmp_file_prefix + '_mds42_refactored.meta')
    with open(metadata_output_file) as metadata_output_fh:
        metadata = pickle.load(metadata_output_fh)

    return (genome_record, metadata)
