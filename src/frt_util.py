"""
Utility methods for inserting FRT sites used for
recombinase-mediated cassette exchange.

This is generally done as the final step before selecting a region of the
genome to send off for synthesis.

http://en.wikipedia.org/wiki/FLP-FRT_recombination
"""

import copy
import os

from Bio import SeqIO

from biopython_util import get_genome_record
from biopython_util import insert_sequence
from refactor_config import OUTPUT_DIR


FRT_WT = 'GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC'

FRT_1 = 'GAAGTTCCTATTCTGTACAAAGTATAGGAACTTC'

FRT_2 = 'GAAGTTCCTATTCGGTACCAAGTATAGGAACTTC'

FRT_3 = 'GAAGTTCCTATTCTCGAGAAAGTATAGGAACTTC'

# Type used for genbank annotation.
FRT_FEATURE_TYPE = 'frt_insert'


def insert_frt_site(
        genome_record,
        upstream_frt_seq,
        upstream_insert_pos,
        downstream_frt_seq,
        downstream_insert_pos,
        feature_id_prefix='',
        upstream_validation_seq=None,
        downstream_validation_seq=None):
    """Inserts FRT sites at the given position.

    Args:
        genome_record: The genome_record to start with. This is not mutated
            in this method.
        upstream_frt: One of FRT_OPTIONS.
        upstream_insert_pos: Insert position for the upstream FRT site.
        downstream_frt: One of FRT_OPTIONS.
        downstream_insert_pos: Insert position for the downstream FRT site.
            NOTE: This is the insert position before the upstream is inserted.
                This method will account for that.
        feature_id_prefix: Prefix for the BioPython feature id.
            Recommended: 'seg_3' for the upstream FRT site of
            segment 3.
        upstream_validation_seq: The next n bases after the start of the
            upstream FRT site. This is used as a sanity check to make sure the
            FRT is being inserted in the right place.
        downstream_validation_seq: The n bases before the downstream site.
            Used for sanity checking the insertion.

    Returns:
        A modified SeqRecord with the change.
    """
    # Maybe check the validation bases.
    if upstream_validation_seq:
        actual_upstream_seq = str(
                genome_record.seq[upstream_insert_pos:
                    upstream_insert_pos + len(upstream_validation_seq)])
        assert upstream_validation_seq == actual_upstream_seq, (
                "Actual %s" % actual_upstream_seq)
    if downstream_validation_seq:
        actual_downstream_seq = str(
                genome_record.seq[
                        downstream_insert_pos - len(downstream_validation_seq):
                        downstream_insert_pos])
        assert downstream_validation_seq == actual_downstream_seq, (
                "Actual %s" % actual_downstream_seq)


    # Adjust the downstream insert position to account for the upstream
    # one being inserted first.
    adjusted_downstream_insert_pos = (
            downstream_insert_pos + len(upstream_frt_seq))

    # Make a copy of the record so we can modify it.
    updated_genome_record = copy.deepcopy(genome_record)

    # Insert the upstream frt.
    updated_genome_record = insert_sequence(
            updated_genome_record,
            upstream_frt_seq,
            upstream_insert_pos,
            safe_features=[],
            insert_feature_type=FRT_FEATURE_TYPE,
            insert_feature_id=feature_id_prefix + '_upstream_frt',
            insert_feature_strand=1)

    # Insert the downstream frt.
    updated_genome_record = insert_sequence(
            updated_genome_record,
            downstream_frt_seq,
            adjusted_downstream_insert_pos,
            safe_features=[],
            insert_feature_type=FRT_FEATURE_TYPE,
            insert_feature_id=feature_id_prefix + '_downstream_frt',
            insert_feature_strand=1)

    return updated_genome_record


if __name__ == '__main__':
    # Test params for seg1 (we did this manually).
    GENOME_RECORD = get_genome_record(
            '../data/completed_segments/seg1/'
            '2013_02_15_04_54_49_seg_1_final.gbk')
    UPSTREAM_FRT = FRT_WT
    UPSTREAM_INSERT_POS = 52749
    UPSTREAM_NEXT_10_BASES = 'GGAGATAGCC'
    FEATURE_ID_PREFIX = 'seg1'
    DOWNSTREAM_FRT = FRT_1[1:]
    DOWNSTREAM_INSERT_POS = 100864
    DOWNSTREAM_PRECEDING_10_BASES = 'AATTTTTATG'

    UPDATED_GENOME_RECORD = insert_frt_site(
            GENOME_RECORD,
            UPSTREAM_FRT,
            UPSTREAM_INSERT_POS,
            DOWNSTREAM_FRT,
            DOWNSTREAM_INSERT_POS,
            FEATURE_ID_PREFIX,
            UPSTREAM_NEXT_10_BASES,
            DOWNSTREAM_PRECEDING_10_BASES)

    # Write the resulting genome.
    GENOME_OUTPUT_FILE = os.path.join(
            OUTPUT_DIR, 'mds42_full_test_frt.gbk')
    with open(GENOME_OUTPUT_FILE, 'w') as output_fh:
        SeqIO.write(UPDATED_GENOME_RECORD, output_fh, 'genbank')

