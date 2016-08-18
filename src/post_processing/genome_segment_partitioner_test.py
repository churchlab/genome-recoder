"""
Tests for preparing part of a genome for synthesis by Gen9 .
"""

import unittest

from biopython_util import get_genome_record
from genome_segment_partitioner import PartitionedSequence


# TODO: Working on public release. Fix these tests.

# class TestPartitionSegment(unittest.TestCase):
#     """Tests for the main partitioning method.
#     """

#     def test_partition_segment__gen9_3kb(self):
#         """Test partitioning with Gen9 params.
#         """
#         GENOME_RECORD = get_genome_record('../data/mds42_full.gbk')
#         START_POSITION = 50000
#         END_POSITION = 100000
#         FRAGMENT_SIZE = 3000
#         MAX_FRAGMENT_SIZE = 3100
#         MAX_WIGGLE = 50

#         forward_seq = str(GENOME_RECORD.seq[START_POSITION:END_POSITION])
#         partitioned_seq_obj = PartitionedSequence(
#                 forward_seq, FRAGMENT_SIZE, MAX_FRAGMENT_SIZE, MAX_WIGGLE)

#         # Test the number of fragments.
#         EXPECTED_FRAGMENTS = 17
#         actual_num_fragments = (
#                 partitioned_seq_obj.get_num_fragments_based_on_seq_length())
#         self.assertEqual(EXPECTED_FRAGMENTS, actual_num_fragments)

#         # Test finding a valid partitioning.
#         fragment_list = partitioned_seq_obj.get_fragment_list()
#         self.assertEqual(EXPECTED_FRAGMENTS, len(fragment_list))


#     def test_partition_segment__agilent_1kb(self):
#         """Test partitioning with Agilent params.
#         """
#         GENOME_RECORD = get_genome_record('../data/mds42_full.gbk')
#         START_POSITION = 50000
#         END_POSITION = 100000
#         FRAGMENT_SIZE = 1000
#         MAX_FRAGMENT_SIZE = 1200
#         MAX_WIGGLE = 100

#         forward_seq = str(GENOME_RECORD.seq[START_POSITION:END_POSITION])
#         partitioned_seq_obj = PartitionedSequence(
#                 forward_seq, FRAGMENT_SIZE, MAX_FRAGMENT_SIZE, MAX_WIGGLE)

#         # Test the number of fragments.
#         EXPECTED_FRAGMENTS = 46
#         actual_num_fragments = (
#                 partitioned_seq_obj.get_num_fragments_based_on_seq_length())
#         self.assertEqual(EXPECTED_FRAGMENTS, actual_num_fragments)

#         # Test finding a valid partitioning.
#         fragment_list = partitioned_seq_obj.get_fragment_list()
#         self.assertEqual(EXPECTED_FRAGMENTS, len(fragment_list))


# if __name__ == '__main__':
#     unittest.main()
