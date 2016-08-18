"""
Utility methods for performing a diff between genomes.
"""

import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from biopython_util import get_feature_by_id
from biopython_util import get_genome_record

COMMENT_CHAR = '#'

def safe_opener(first, second, file_format):
    """Returns tuple (first_seq_record, second_seq_record).

    Checks that they are the same length.
    """
    if isinstance(first, SeqRecord) and isinstance(second, SeqRecord):
        return (first, second)

    if file_format == 'genbank':
        print 'Reading %s ...' % os.path.split(first)[1]
        first_genome_record = get_genome_record(first,
                use_old_id_strategy=True)
        print 'Reading %s ...' % os.path.split(second)[1]
        second_genome_record = get_genome_record(second,
                use_old_id_strategy=True)
    elif file_format == 'fasta':
        with open(first) as first_fh:
            first_genome_record = SeqIO.read(first_fh, file_format)

        with open(second) as second_fh:
            second_genome_record = SeqIO.read(second_fh, file_format)
    else:
        raise AssertionError("File format %s not supported." % file_format)

    # For now we require them to be the same size as we are not setup to handle
    # insertions or deletions.
    assert len(first_genome_record) == len(second_genome_record)

    return (first_genome_record, second_genome_record)


def find_mismatches_between_same_size_genomes(first, second, report_file,
        file_format='genbank'):
    """Reports all the positions that don't match between
    the first and second genome.
    """
    (first_genome_record, second_genome_record) = safe_opener(first, second,
            file_format)

    first_genome_record_seq_str = str(first_genome_record.seq)
    second_genome_record_seq_str = str(second_genome_record.seq)

    with open(report_file, 'w') as report_fh:
        for base_index in range(len(first_genome_record_seq_str)):
            if (not first_genome_record_seq_str[base_index] ==
                    second_genome_record_seq_str[base_index]):
                report_fh.write('MISMATCH at %d. First: %s. Second: %s\n' % (
                        base_index, first_genome_record_seq_str[base_index],
                        second_genome_record_seq_str[base_index]))


def find_features_differing_in_position(first, second, report_file,
        file_format='genbank', ignore_ids_source=None):
    """Identify features that are different in position, or altogether missing
    from one genome but not the other.
    """
    (first_genome_record, second_genome_record) = safe_opener(first, second,
            file_format)

    # Grab the features we are ignoring.
    ignore_ids = set([])
    if ignore_ids_source is not None:
        with open(ignore_ids_source) as ignore_ids_fh:
            for line in ignore_ids_fh:
                clean_line = line.strip()
                if len(clean_line) > 0 and clean_line[0] != COMMENT_CHAR:
                    ignore_ids.add(line.strip())

    print 'Starting analysis ...'
    with open(report_file, 'w') as report_fh:

        # First, figure out feature ids that are different btw the two.
        first_genome_feature_id_set = set([feature.id for feature in
                first_genome_record.features])
        second_genome_feature_id_set = set([feature.id for feature in
                second_genome_record.features])
        in_first_but_not_second = (first_genome_feature_id_set -
                second_genome_feature_id_set - ignore_ids)
        in_second_but_not_first = (second_genome_feature_id_set -
                first_genome_feature_id_set - ignore_ids)

        report_fh.write('In %s but missing from %s:\n' % (
                first_genome_record.name, second_genome_record.name))
        for feature_id in in_first_but_not_second:
            report_fh.write(feature_id + '\n')

        report_fh.write('\n\n')

        report_fh.write('In %s but missing from %s:\n' % (
                second_genome_record.name, first_genome_record.name))
        for feature_id in in_second_but_not_first:
            report_fh.write(feature_id + '\n')


if __name__ == '__main__':
    find_features_differing_in_position(
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/rec1_c321d.genbank',
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/C321(F1,M1,Ao)_I3_not_strict.gb',
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/diff_report_2.txt',
            'genbank',
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/diff_safe_to_ignore.txt')
