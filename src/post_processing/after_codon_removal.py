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
More actions to take after removing all forbidden codons.

These steps are in development so they are done more-or-less manually
in this script. Steps include:
    * remove homopolymer runs
    * remove restriction sites
    * add end pieces (e.g. FRT sites)

The only remaining step after this is to either perform Gen9 partitioning
or just excise the segment for synthesis by SGI. The methods to do this are
not called in this module.
"""

import ast
import copy
import string

import pandas as pd

from biopython_util import swap_unique_seq
from frt_util import insert_frt_site
from homopolymers import remove_homopolymer_runs
from gc_content_fixer import fix_gc_content
from gc_content_fixer import GCContentConstraints
from restriction_sites import remove_restriction_sites
from refactor_checker import check_recoding_is_complete
from refactor_config import ORIGINAL_GENOME_RECORD
from refactor_config import PRETTY_PRINTER
from refactor_config import RESTRICTION_ENZYME_SITES_TO_REMOVE
from refactor_context import RefactorContext


def perform_final_steps(refactor_context, seg_start, seg_end,
        upstream_flanking_seq=None, downstream_flanking_seq=None,
        validation_start_seq=None, validation_end_seq=None,
        ignore_problems_in_feature_ids=[], report_prefix=None):
    """Updates the contained genome_record after mutating it, including:
        * remove homopolymer runs
        * remove restriction sites
        * (optional) add end pieces (e.g. FRT sites)

    Args:
        refactor_context: The RefactorContext.
        seg_start: The first position (pythonic) in the genome_record contained
            within refactor_context for the segment.
        seg_end: End position (pythonic) for the segment.
        upstream_flanking_seq: Sequence to insert at the head of the segment.
        downstream_flanking_seq: Sequence to insert at the tail of the segment.
        validation_start_seq: Optional sequence at the start of the segment
            to sanity check the start position. We lack a ui.
        validation_end_seq: Optional sequence at the end of the segment
            to sanity check the end position.
        ignore_problems_in_feature_ids: Feature ids that the client
            is aware may have problems so that we can ignore.

    Returns:
        An updated SeqRecord reflecting changes.
    """
    updated_genome_record = copy.deepcopy(refactor_context.get_genome_record())

    # Check features are conserved before we start.
    check_recoding_is_complete(ORIGINAL_GENOME_RECORD, updated_genome_record,
            ignore_problems_in_feature_ids=ignore_problems_in_feature_ids,
            interval=(seg_start, seg_end))

    orig_seq = str(updated_genome_record.seq)

    if validation_start_seq:
        assert validation_start_seq == orig_seq[
                seg_start:seg_start + len(validation_start_seq)]

    if validation_end_seq:
        assert validation_end_seq == orig_seq[
                seg_end - len(validation_end_seq):seg_end]

    updated_refactor_context = RefactorContext(updated_genome_record)

    # Fix GC content.
    GC_CONTENT_CONSTRAINT_OBJ = GCContentConstraints()
    updated_genome_record = fix_gc_content(
            refactor_context,
            GC_CONTENT_CONSTRAINT_OBJ,
            start_bound=seg_start,
            end_bound=seg_end,
            debug=False)
    updated_refactor_context.set_genome_record(updated_genome_record)

    # Remove homopolymer runs.
    remove_homopolymer_result = remove_homopolymer_runs(
            updated_refactor_context, start_bound=seg_start, end_bound=seg_end,
            report_prefix=report_prefix)
    updated_genome_record = remove_homopolymer_result['updated_genome_record']
    flagged_h_runs = remove_homopolymer_result['flagged']
    updated_refactor_context.set_genome_record(updated_genome_record)
    print 'Flagged homopolyer runs:'
    PRETTY_PRINTER.pprint(flagged_h_runs)

    # Remove restriction sites.
    remove_res_sites_result = remove_restriction_sites(
            updated_refactor_context, RESTRICTION_ENZYME_SITES_TO_REMOVE,
            start_bound=seg_start, end_bound=seg_end,
            report_prefix=report_prefix)
    updated_genome_record = remove_res_sites_result['updated_genome_record']
    updated_refactor_context.set_genome_record(updated_genome_record)
    flagged_res_sites = remove_res_sites_result['flagged']
    print 'Flagged restriction sites:'
    PRETTY_PRINTER.pprint(flagged_res_sites)

    # Generate the GC content report after all fixes are done.
    gc_report_file = report_prefix + 'gc_content.csv'
    updated_genome_record = fix_gc_content(
            updated_refactor_context,
            GC_CONTENT_CONSTRAINT_OBJ,
            start_bound=seg_start,
            end_bound=seg_end,
            debug=True,
            report_file=gc_report_file)

    # Check features are conserved.
    print 'Checking translation/rna/forbidden codons ...'
    check_recoding_is_complete(ORIGINAL_GENOME_RECORD, updated_genome_record,
            ignore_problems_in_feature_ids=ignore_problems_in_feature_ids,
            interval=(seg_start, seg_end))

    # Maybe insert FRT sites.
    if upstream_flanking_seq or downstream_flanking_seq:
        updated_genome_record = insert_frt_site(
                updated_genome_record,
                upstream_flanking_seq,
                seg_start,
                downstream_flanking_seq,
                seg_end,
                feature_id_prefix='seg2',
                upstream_validation_seq=validation_start_seq,
                downstream_validation_seq=validation_end_seq)

    return updated_genome_record


def handle_fixes_in_file(genome_record, fixes_file):
    """Handles fixes in provided csv file.

    User should provide a file with origina/new sequence of ~30 bp that is
    unique to within some epsilon, default 100 bp. See swap_unique_seq() params.

    Args:
        genome_record: Mutable genome_record.
        fixes_file: Csv with at least the following columns:
            * original: Original sequence
            * new: New version of the sequence
            * interval: Approximate interval of the target of the fix.
    """
    print "handling fixes in %s" % fixes_file
    fixes_df = pd.read_csv(fixes_file)
    assert 'original' in fixes_df
    assert 'new' in fixes_df
    assert 'interval' in fixes_df
    fixes_df['interval'] = fixes_df['interval'].apply(ast.literal_eval)

    for idx, row in fixes_df.iterrows():
        if isinstance(row['original'], float):
            # NaN
            continue
        try:
            swap_unique_seq(
                    genome_record,
                    string.upper(row['original']),
                    string.upper(row['new']),
                    fuzzy_region=row['interval'])
        except AssertionError:
            print 'WARNING: Could not find %s' % string.upper(row['original'])
            # Probably couldn't find it because overlaps, just continue on.
            continue


if __name__ == '__main__':
    from Bio import SeqIO

    from biopython_util import get_genome_record
    from refactor_context import RefactorContext

    genome_record = get_genome_record(
            '../data/completed_segments/seg2/2013_03_06_20_16_04_mds42_refactored.gbk')
    genome_record.name =  genome_record.name[:-3] + 'seg2'
    refactor_context = RefactorContext(genome_record)

    SEG_START = 100863
    SEG_END = 148475
    UPSTREAM_FLANKING_SEQ = 'CAGCCTTGTTTCGCCAGAATGCCAGTCAGCATAAGGGAGAGCTCAAGGCAGAAGTTCCTATTCCGAAGTTCCTATTCTCATATAAGTATAGGAACTTC'
    DOWNSTREAM_FLANKING_SEQ = 'CCTGTTGACAATTAATCATCGGCATAGTATATCGGCATAGTATAATACGACAAGGTGAGGAACTAAACCCAGGAGGCAGATCATGAGTCTGAAAGAAAAAACACAATCTCTGTTTGCCAACGCATTTGGCTACCCTGCCACTCACACCATTCAGGCGCCTGGCCGCGTGAATTTGATTGGTGAACACACCGACTACAACGACGGTTTCGTTCTGCCCTGCGCGATTGATTATCAAACCGTGATCAGTTGTGCACCACGCGATGACCGTAAAGTTCGCGTGATGGCAGCCGATTATGAAAATCAGCTCGACGAGTTTTCCCTCGATGCGCCCATTGTCGCACATGAAAACTATCAATGGGCTAACTACGTTCGTGGCGTGGTGAAACATCTGCAACTGCGTAACAACAGCTTCGGCGGCGTGGACATGGTGATCAGCGGCAATGTGCCGCAGGGTGCCGGGTTAAGTTCTTCCGCTTCACTGGAAGTCGCGGTCGGAACCGTATTGCAGCAGCTTTATCATCTGCCGCTGGACGGCGCACAAATCGCGCTTAACGGTCAGGAAGCAGAAAACCAGTTTGTAGGCTGTAACTGCGGGATCATGGATCAGCTAATTTCCGCGCTCGGCAAGAAAGATCATGCCTTGCTGATCGATTGCCGCTCACTGGGGACCAAAGCAGTTTCCATGCCCAAAGGTGTGGCTGTCGTCATCATCAACAGTAACTTCAAACGTACCCTGGTTGGCAGCGAATACAACACCCGTCGTGAACAGTGCGAAACCGGTGCGCGTTTCTTCCAGCAGCCAGCCCTGCGTGATGTCACCATTGAAGAGTTCAACGCTGTTGCGCATGAACTGGACCCGATCGTGGCAAAACGCGTGCGTCATATACTGACTGAAAACGCCCGCACCGTTGAAGCTGCCAGCGCGCTGGAGCAAGGCGACCTGAAACGTATGGGCGAGTTGATGGCGGAGTCTCATGCCTCTATGCGCGATGATTTCGAAATCACCGTGCCGCAAATTGACACTCTGGTAGAAATCGTCAAAGCTGTGATTGGCGACAAAGGTGGCGTACGCATGACCGGCGGCGGATTTGGCGGCTGTATCGTCGCGCTGATCCCGGAAGAGCTGGTGCCTGCCGTACAGCAAGCTGTCGCTGAACAATATGAAGCAAAAACAGGTATTAAAGAGACTTTTTACGTTTGTAAACCATCACAAGGAGCAGGACAGTGCTGAAAAAAAAAACCCCGCCCCTGACAGGGCGGGGTTTTTTTTGAAGTTCCTATTCCGAAGTTCCTATTCTATCAGAAGTATAGGAACTTCAGTGCGGATTTCGTATTTGCAGCTCGTCAGTACTTTCAGAATCATGGCCT'
    VALIDATION_START_SEQ = 'GAGGCCGACGATGATTACGGCCTCAGG'
    VALIDATION_END_SEQ = 'TTAATCATTTGACGTCCCTTGT'

    updated_genome_record = perform_final_steps(refactor_context, SEG_START,
            SEG_END, UPSTREAM_FLANKING_SEQ, DOWNSTREAM_FLANKING_SEQ,
            VALIDATION_START_SEQ, VALIDATION_END_SEQ)

    genome_output_file = (
            '../data/completed_segments/seg2/seg2_final_with_flanking.gbk')
    with open(genome_output_file, 'w') as output_fh:
        SeqIO.write(updated_genome_record, output_fh, 'genbank')
