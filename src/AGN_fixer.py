"""
Gonna fix some AGNs.
"""

import csv
import os
import re

from Bio import SeqIO

from biopython_util import get_feature_by_id
from biopython_util import get_genome_record
from biopython_util import swap_feature_codon_at_position
from paths import DATA_DIR
from refactor_config import ORIGINAL_CODON_USAGE_MEMEX
from refactor_config import ORIGINAL_GENOME_RECORD

MDS42_RECORD = ORIGINAL_GENOME_RECORD

SEG5_2MB_DIR = os.path.join(DATA_DIR, 'completed_segments', 'seg5-seg48')
RECODED_PATH = os.path.join(SEG5_2MB_DIR,
        '2013_08_30.seg5-seg48.motifs_removed.genbank')

OUTFILE = os.path.join(SEG5_2MB_DIR,
        '2013_08_30.seg5-seg48.motifs_removed.AGN_fixed.genbank')

AGN_DEBUG_FILE = os.path.join(SEG5_2MB_DIR, 'AGN_debug_FINAL_MJL.csv')

def main():
    record = get_genome_record(RECODED_PATH)

    with open(AGN_DEBUG_FILE) as agn_debug_fh:
        reader = csv.DictReader(agn_debug_fh)
        for row in reader:
            new_codon = row['Codon']
            if len(new_codon) != 3:
                continue

            feature_id = row['ID']
            print feature_id

            # Identify codon positions from original record.
            sequence = row['Sequence']
            if len(sequence) < 14:
                # Too hard.
                continue
            orig_codon = sequence[8:11]
            if not re.match(r'[A-Z]{3}', orig_codon):
                continue

            # Identify which codon index the position is.
            seq_ending_with_codon = sequence[0:3].upper() + sequence[4:7].upper() + orig_codon
            print seq_ending_with_codon

            orig_feature = get_feature_by_id(MDS42_RECORD, feature_id)
            orig_feature_seq = str(orig_feature.extract(MDS42_RECORD.seq))

            new_feature = get_feature_by_id(record, feature_id)
            new_feature_seq = str(new_feature.extract(record.seq))

            # Find the last occurrence.
            last_pos = None
            for match in re.finditer(seq_ending_with_codon, orig_feature_seq):
                last_pos = match.start()
            print last_pos
            codon_pos = last_pos + 6
            print 'expect', orig_feature_seq[codon_pos:codon_pos + 3]

            # Translate to codon index.
            codon_index = codon_pos / 3

            # Make sure codon at that position in recoded record is synonymous.
            codon_at_index = new_feature_seq[codon_pos:codon_pos + 3]
            print codon_at_index
            syn = ORIGINAL_CODON_USAGE_MEMEX.get_synonymous_codons(codon_at_index)
            print syn
            assert orig_codon in syn

            # Make the change in the target record.
            swap_feature_codon_at_position(
                    record, new_feature.id, codon_pos, codon_at_index, new_codon)

            # Maybe do the next one.
            maybe_next = row['next'].upper()
            if len(maybe_next) == 3:
                syn = ORIGINAL_CODON_USAGE_MEMEX.get_synonymous_codons(maybe_next)
                codon_at_index_next = new_feature_seq[codon_pos + 3:codon_pos + 6]
                assert codon_at_index_next in syn
                swap_feature_codon_at_position(
                        record, new_feature.id, codon_pos + 3, codon_at_index_next, maybe_next)

    with open(OUTFILE, 'w') as fh:
        SeqIO.write(record, fh, 'genbank')


if __name__ == '__main__':
    main()
