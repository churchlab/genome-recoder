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
Utilities for working with the CD-Hit sequence aligner.

See: http://weizhong-lab.ucsd.edu/cd-hit/
"""

import os
import subprocess

from cogent.app import cd_hit as cogent_cd_hit_parser

from refactor_config import TMP_DATA_DIR


DEFAULT_THRESHOLD = 0.6

DEFAULT_WORD_SIZE = 4


class CDHitWrapper(object):
    """Wrapper object for the CD-Hit sequencer aligner.
    """

    def __init__(self, debug=False, **kwargs):
        """Constructor.

        Args:
            debug: Whether to run in debug mode (more noise output).
            kwargs: Optional config arguments. Currently supported:
                * threshold
                * word_size (acceptable values depend on threshold).
        """
        self.debug = debug
        self.threshold = kwargs.get('threshold', DEFAULT_THRESHOLD)
        self.word_size = kwargs.get('word_size', DEFAULT_WORD_SIZE)


    def detect_high_homology_pairs(self, seq_list):
        """Given a list of sequences, performs an all-to-all alignment
        and returns pairs of sequences (as indexes in the original list),
        which have a higher than desired homology.

        Args:
            seq_list: List of sequence strings. Order matters in that the
                pairs of indeces returned are relative to the ordering of the
                list.

        Returns:
            A List<Set<int>> of indeces that corresponding to sequences that
            were found to have homology above the given threshold with each
            other.
        """
        # Make sure the tmp directory
        if not os.path.exists(TMP_DATA_DIR):
            os.mkdir(TMP_DATA_DIR)

        input_filename = os.path.join(TMP_DATA_DIR, 'test_cdhit_input.fasta')
        output_prefix = os.path.join(TMP_DATA_DIR, 'test_cdhit_output')

        # Write the list of sequences to a file in .fasta format, as expected
        # by the cd-hit command.
        with open(input_filename, 'w') as input_fh:
            for seq_index in range(len(seq_list)):
                # For each sequence write the following two lines.
                # >0
                # ATGAGATAGTA
                seq = seq_list[seq_index]
                input_fh.write('>' + str(seq_index) + '\n' + str(seq) + '\n')

        # Call cd-hit.
        cmd = [
            'cd-hit',
            '-i', input_filename,
            '-o', output_prefix,
            '-c', str(self.threshold),
            '-n', str(self.word_size),
        ]
        if self.debug:
            stdout = None
        else:
            devnull = open(os.devnull, 'wb')
            stdout = devnull
        subprocess.call(cmd, stdout=stdout)

        # Parse the results.
        cluster_output_filename = output_prefix + '.clstr'
        with open(cluster_output_filename) as cluster_output_fh:
            cluster_lines = [line.rstrip() for line in cluster_output_fh]

        cluster_list = cogent_cd_hit_parser.parse_cdhit_clstr_file(
                cluster_lines)

        # The homology conflicts are clusters with more than one element.
        homology_clusters = []
        for cluster in cluster_list:
            if len(cluster) > 1:
                homology_clusters.append(set(
                    [int(seq_index_str) for seq_index_str in cluster]))

        return homology_clusters
