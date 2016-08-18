"""
Muddle up all the RBS cp regions.
"""

import csv
import os
import re

from Bio import SeqIO

from biopython_util import get_feature_by_id
from biopython_util import get_genome_record
from biopython_util import InsertType
from biopython_util import update_seq_record_feature
from codon_replacer_util import replace_codons_in_single_feature
from refactor_context import RefactorContext

from paths import DATA_DIR

SEG5_2MB_DIR = os.path.join(DATA_DIR, 'completed_segments', 'seg5-seg48')
RECODED_PATH = os.path.join(SEG5_2MB_DIR,
        '2013_08_30.seg5-seg48.motifs_removed.AGN_fixed.genbank')

ID_ROOT = 'ECMDS42_'

OUTFILE = os.path.join(SEG5_2MB_DIR,
        '2013_08_30.seg5-seg48.motifs_removed.AGN_fixed.muddled.genbank')

AGN_DEBUG_FILE = os.path.join(SEG5_2MB_DIR, 'AGN_debug_FINAL_MJL.csv')

def main():
    source_ids_to_muddle = []
    with open(AGN_DEBUG_FILE) as agn_debug_fh:
        reader = csv.DictReader(agn_debug_fh)
        for row in reader:
            if row['Separate'] == '1':
                source_ids_to_muddle.append(row['ID'])
    print source_ids_to_muddle

    record = get_genome_record(RECODED_PATH)
    refactor_context = RefactorContext(record)

    for source_feature_id in source_ids_to_muddle:
        source_feature = get_feature_by_id(record, source_feature_id)
        muddle_end(source_feature, record, refactor_context, 20)

    # rbs_cp_features = [feature for feature in record.features if
    #         feature.type == InsertType.FIX_OVERLAP_RBS_COPY]

    # overlap_head_cp_features = [feature for feature in record.features if
    #         feature.type == InsertType.FIX_OVERLAP_HEAD_COPY]

    # print 'rbs', len(rbs_cp_features)
    # print 'head', len(overlap_head_cp_features)

    # head_cp_feature_ids = set()
    # for head_feature in overlap_head_cp_features:
    #     source_feature_id = re.match(r'(?P<feature_id>.*)_' + InsertType.FIX_OVERLAP_HEAD_COPY, head_feature.id).group('feature_id')
    #     head_cp_feature_ids.add(source_feature_id)

    # count = 0
    # for rbs_cp_feature in rbs_cp_features:
    #     downstream = False
    #     match = re.match(r'(?P<feature_id>.*)_upstream_' + InsertType.FIX_OVERLAP_RBS_COPY, rbs_cp_feature.id)
    #     if not match:
    #         downstream = True
    #         match = re.match(r'(?P<feature_id>.*)_downstream_' + InsertType.FIX_OVERLAP_RBS_COPY, rbs_cp_feature.id)
    #     source_feature_id = match.group('feature_id')
    #     if source_feature_id in head_cp_feature_ids:
    #         print 'HAS HEAD_CP', source_feature_id
    #         continue
    #     source_feature = get_feature_by_id(record, source_feature_id)
    #     num_part = re.match(r'.*_(?P<num>[0-9]+)', source_feature_id).group('num')

    #     if source_feature.strand == 1:
    #         assert source_feature.location.start > rbs_cp_feature.location.start
    #         actual_source_id = ID_ROOT + str(int(num_part) - 1)
    #     else:
    #         assert source_feature.location.start < rbs_cp_feature.location.start
    #         actual_source_id = ID_ROOT + str(int(num_part) + 1)
    #     try:
    #         actual_source_feature = get_feature_by_id(record, actual_source_id)
    #     except:
    #         print 'NOT MUDDLING', rbs_cp_feature.id
    #         continue

    #     if actual_source_feature.strand != rbs_cp_feature.strand:
    #         print 'NOT MUDDLING', rbs_cp_feature.id
    #         continue

    #     print 'MUDDLING', rbs_cp_feature.id
    #     muddle_end(actual_source_feature, record, refactor_context, len(rbs_cp_feature))

    #     count += 1

    # print 'COUNT', count

    with open(OUTFILE, 'w') as fh:
        SeqIO.write(record, fh, 'genbank')


def muddle_end(source_feature, genome_record, refactor_context, size_to_modify):
    num_codons = len(source_feature) / 3

    first_codon_to_modify = num_codons - size_to_modify / 3
    source_feature_seq = source_feature.extract(genome_record.seq)

    avoid_codons_in_positions = {}
    for codon_index in range(first_codon_to_modify, num_codons):
        codon = str(source_feature_seq[codon_index * 3 : codon_index * 3 + 3])
        avoid_codons_in_positions[codon_index] = codon

    # Perform the fix.
    result = replace_codons_in_single_feature(
            refactor_context,
            source_feature.id,
            start_codon_index=first_codon_to_modify,
            avoid_codons_in_positions=avoid_codons_in_positions)
    assert str(result['orig_feature_seq']) != str(result['new_feature_seq'])
    assert result['is_success'], "Resolving homology not successful."

    update_seq_record_feature(
            genome_record,
            source_feature.id,
            result
    )


if __name__ == '__main__':
    main()
