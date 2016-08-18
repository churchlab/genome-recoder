"""
Module with specific alignment methods that we use.

Also tests to calibrate the parameters that we want.

One place we use pairwise aligners is in checking for homology among
overhangs when partitioning a ~48kb sequence into fragments to send off to
Gen9, which requires pieces of size 3kb (and smaller).
"""

from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment


global MATCH
global MISS
global GAP_START
global GAP_INC

# Homology params
MATCH = 2
MISS = -1
GAP_START = -1
GAP_INC = -1


def calc_best_pairwise_local_align_score(seq1, seq2, score_only=True):
    """Find the highest score from a local alignment between the two sequences.

    For now we are using the Bioython pairwise2 module.

    TODO: Make sense of the magnitude of the score.
    """
    return align.localms(seq1, seq2,
            MATCH, MISS, GAP_START, GAP_INC, score_only=score_only)


def test_bio_pairwise2():
    # Sequences uesd for testing.
    simple_seq = 'AAAAAAAAAA' + 'AAAAAAAAAA' + 'AAAAAAAAAA' + 'AAAAAAAAAA' + 'AAAAAAAAAA'
    # Align a repeating sequence to itself.
    print calc_best_pairwise_local_align_score(simple_seq, simple_seq)

    # Try aligning with a mismatch mid-run.
    simple_seq_2 = 'AAAAAAAAAA' + 'TTTTTTTTTT' + 'AAAAAAAAAA' + 'AAAAAAAAAA' + 'AAAAAAAAAA'
    print calc_best_pairwise_local_align_score(simple_seq, simple_seq_2)

    # Try completely non-homologous sequence.
    simple_seq_3 = 'CCCCCCCCCC' + 'CCCCCCCCCC' + 'CCCCCCCCCC' + 'CCCCCCCCCC' + 'CCCCCCCCCC'
    print calc_best_pairwise_local_align_score(simple_seq, simple_seq_3)

    # Partial alignment seems to start.
    simple_seq_4 = 'AAAAAAAAAA' + 'CCCCCCCCCC' + 'CCCCCCCCCC' + 'CCCCCCCCCC' + 'CCCCCCCCCC'
    print calc_best_pairwise_local_align_score(simple_seq, simple_seq_4)

    # Partial alignment seems to start.
    simple_seq_5 = 'AAAACCCAAA' + 'CCCCCCCCCC' + 'CCCCCCCCCC' + 'AAAAAAAAAA' + 'CCCCCCCCCC'
    print calc_best_pairwise_local_align_score(simple_seq, simple_seq_5)

    # Flip-flop almost match.
    simple_seq_6 = 'ATATATATAT' + 'ATATATATAT' + 'ATATATATAT' + 'ATATATATAT' + 'ATATATATAT'
    print calc_best_pairwise_local_align_score(simple_seq, simple_seq_6)



if __name__ == '__main__':
    from Bio import SeqIO
    overhang_seqs = []
    with open('tmp/test_cdhit_input.fasta') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            overhang_seqs.append(str(record.seq))


    analysis_file = 'tmp/gen9_fragment_homology_analysis_combo.txt'

    with open(analysis_file, 'w') as fh:
        fh.write('Homology sanity check using python pairwise2\n\n')

        for gap_iter in [-1, -2]:
            for mismatch_iter in [-1, -2, -3, -4]:
                GAP_START = gap_iter
                GAP_INC = gap_iter
                MISS = mismatch_iter

                # Map score to pair that produced the score so we can print out the
                # strongest homologies and visually inspect them.
                score_to_pair_map = {}
                for i in range(len(overhang_seqs)):
                    for j in range(i + 1, len(overhang_seqs)):
                        seq1 = overhang_seqs[i]
                        seq2 = overhang_seqs[j]
                        score = calc_best_pairwise_local_align_score(seq1, seq2)
                        score_to_pair_map[score] = (i, j)
                top_10_scores = sorted(score_to_pair_map.keys(), reverse=True)[:10]

                fh.write('\n========================================'
                        '========================================\n'
                        '========================================'
                        '========================================\n\n')

                fh.write('Params:\n')
                fh.write('\tMATCH: %d\n' % MATCH)
                fh.write('\tMISS: %d\n' % MISS)
                fh.write('\tGAP_START: %d\n' % GAP_START)
                fh.write('\tGAP_INC: %d\n' % GAP_INC)

                # For the worst 10, printout a few examples of bad homology.
                for score in top_10_scores:
                    (i, j) = score_to_pair_map[score]
                    seq1 = overhang_seqs[i]
                    seq2 = overhang_seqs[j]
                    fh.write('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
                            '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n')
                    fh.write('SEQ_1: ' + seq1 + '\n')
                    fh.write('SEQ_2: ' + seq2 + '\n\n')
                    alignments = calc_best_pairwise_local_align_score(
                            seq1, seq2, score_only=False)
                    a = alignments[0]
                    fh.write(format_alignment(*a))
