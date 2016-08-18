"""
Logic for partitioning a genome for synthesis according to boundaries defined
by operons.
"""

import csv
import re

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation

from biopython_util import add_feature_to_seq_record
from biopython_util import build_gene_synonym_map
from biopython_util import build_gene_to_CDS_feature_map
from biopython_util import build_linked_cds_feature_map
from biopython_util import get_feature_gene
from biopython_util import get_feature_by_gene_name
from biopython_util import InsertType
from paths import MG1655_SOURCE
from regulondb import identify_mg1655_partition_positions


class MG1655DescendantPartitionHelper(object):
    """Object that manages translating intervals between MG1655 and descendant
    strains.
    """

    def __init__(self, recoded_genome_record, partition_size_lower_bound,
            partition_size_upper_bound):
        """Constructor.
        """
        self.mg1655_genome_record = SeqIO.read(MG1655_SOURCE, 'genbank')

        self.recoded_genome_record = recoded_genome_record

        # Build a map that allows marching through recoded CDS features.
        self.recoded_linked_cds_feature_map = build_linked_cds_feature_map(
                self.recoded_genome_record)

        # Build a map from all possible genes to synonyms of that gene for
        # MG1655.
        self.mg1655_gene_synonym_map = build_gene_synonym_map(
                self.mg1655_genome_record)

        # Build a map from gene to the next gene that follows spatially.
        # Also cover synonyms among keys.
        self.mg1655_linked_cds_map = build_linked_cds_feature_map(
                self.mg1655_genome_record, self.mg1655_gene_synonym_map)

        # We need a map from intervals to corresponding next intervals
        # and genes.
        self.mg1655_interval_to_gap_interval_map = (
                identify_mg1655_partition_positions(
                        partition_size_lower_bound=partition_size_lower_bound,
                        partition_size_upper_bound=partition_size_upper_bound))

        # Build a map from gene to the interval that immediately follows it.
        self.mg1655_gene_to_preceeding_interval_map = {}
        interval_map_iter_items = (
                self.mg1655_interval_to_gap_interval_map.iteritems())
        for interval, interval_dict in interval_map_iter_items:
            self.mg1655_gene_to_preceeding_interval_map[
                    interval_dict['genes'][0][0]] = interval


    def get_interval_before_gene(self, gene):
        """Returns the interval immediately before the given gene, or
        None if such an interval is not known.
        """
        return self.mg1655_gene_to_preceeding_interval_map.get(gene, None)


    def get_next_interval_options(self, interval):
        """Returns a list of next interval options for the given interval.
        """
        return self.mg1655_interval_to_gap_interval_map[interval][
                'following_intervals']


    def get_gene_following_interval(self, interval):
        """Returns the name of the gene following the given interval.
        """
        return self.mg1655_interval_to_gap_interval_map[interval][
                'genes'][0][0]


    def get_next_mg1655_gene(self, gene):
        """Returns the gene that follows spatially after the current gene.
        """
        return self.mg1655_linked_cds_map.get(gene, None)


    def get_gene_synonyms(self, gene):
        """Returns list of gene synonyms. May be an empty list.
        """
        return self.mg1655_gene_synonym_map.get(gene, [])


def is_insertion_element(gene_name):
    """Light-weight, non-thorough check of whether the gene is an insertion
    element.
    """
    return bool(re.match(r'ins.*', gene_name))


def find_segment_partitions(genome_record, num_segments=None, first_gene=None,
        last_gene=None,
        start_seg_numbering=0, seg_annotation_prefix='seg',
        annotated_genome_record_outfile=None,
        output_segment_report=None,
        partition_size_lower_bound=44000, partition_size_upper_bound=52000):
    """Identify intervals to partition the genome.

    Args:
        genome_record: SeqRecord that is mutated.
        num_segments: Number of segments to carve out.
        first_gene: Name of first gene to start at. We use names rather than
            numbers since the partitioning process is inherently fuzzy.
        last_gene: If provided, last segment will have this as last gene.
        start_seg_numbering: Number for first segment.
        seg_annotation_prefix: Name for the segments. (e.g. seg3).
        annotated_genome_record_file: Optional path to write the annotated
            Genbank out to.
    """
    assert first_gene, "Starting from beginning not implemented yet."

    partition_helper = MG1655DescendantPartitionHelper(genome_record,
            partition_size_lower_bound, partition_size_upper_bound)

    gene_to_feature_map = build_gene_to_CDS_feature_map(genome_record)

    # We use the strategy of using genes as bounds for partioning, rather
    # than exact positions. This will help us communicate when using the two
    # different frames where one is MG1655 and the other is the modified MDS42.
    # We might be able to get away with not having to mess around with liftover.


    TEMP_FIRST_GENE_ERROR_MSG = ("Doing our magic gene name de-duping not " +
            "yet implemented for the first gene.")

    # Initialize the starting gene for the next segment, and the interval
    # immediately preceding the gene.
    first_gene_in_segment = first_gene
    assert first_gene_in_segment in gene_to_feature_map, (
            TEMP_FIRST_GENE_ERROR_MSG)
    mg1655_upstream_interval = partition_helper.get_interval_before_gene(
            first_gene_in_segment)
    assert mg1655_upstream_interval, TEMP_FIRST_GENE_ERROR_MSG
    first_gene_in_segment_feature = gene_to_feature_map[first_gene_in_segment]

    # Compute the initial segment start.
    upstream_interval_size = (mg1655_upstream_interval[1] -
            mg1655_upstream_interval[0])
    recoded_seg_start = (
            first_gene_in_segment_feature.location.start -
            upstream_interval_size / 2)

    # Keep track of segments to write to report.
    segment_report_list = []

    for seg_index in range(num_segments):
        # The general algorithm for each segment:
        #     0. The loop starts with the upstream interval bound provided and
        #        the first gene in the segment, which has been checked to exist
        #        both in MG1655 and MDS42.
        #     1. Using the position of the gene in MDS42, we find the last gene
        #        to be included in the segment, so that segment is at least as
        #        big as the lower bound provided in the arguments.
        #     2. Determine the positions to cut within the intervals.
        #         a. Use mid-point as first approxmation.
        #             i. Eventually need to consider promoter sites.
        #     3. Create annotation for this segment.
        #     3. Prepare for the next iteration of the loop.

        # Exceptions:
        #     * If we get to the end of the genome, then just truncate the
        #       final segment. This naive strategy may end up with a
        #       significantly truncated final segment, so will debug from
        #       there.

        ### 0. Flags used at various steps below.

        # Flag to indicate whether we hit the end of the genome.
        hit_end_of_genome = False

        # Flag to indicate we're synthesizing the last segment. This is used
        # in combination with the last_gene param to allow clients to manually
        # extend the final segment being cut.
        is_final_segment = (seg_index + 1 == num_segments)

        # Whether or not take into account last gene / is_final_segment
        # in this iteration.
        should_heed_last_gene = (bool(last_gene) and is_final_segment)

        ### 1. Find the next gene that is a segment distance away.

        def _check_if_end_of_genome(last_gene_in_segment):
            """Helper function used below that checks whether gene is at end of
            genome.

            Returns:
                Boolean indicating whether gene is at end of genome.

            Raises:
                AssertionError if our heuristic for checking end of genome
                is flawed.
            """
            next_segment_feature = (
                    partition_helper.recoded_linked_cds_feature_map.get(
                            last_gene_in_segment))
            if not next_segment_feature:
                # Should be end of the genome unless there's a gap so assert
                # that last segment size is acceptable and use what we have.
                full_seg_size = (len(genome_record) -
                        first_gene_in_segment_feature.location.start)
                assert full_seg_size < partition_size_upper_bound, (
                        "Can't find next segment feature. Not end of genome "
                        "because last segment would be too big: %d" % (
                                full_seg_size,))
                return True
            return False

        # First we find a start last segment that is the minimal distance
        # away from the start gene to make the segment cutoff. Use last_gene
        # if provided.
        # If we reach the end of the genome in this step, just use the last
        # gene found.
        if should_heed_last_gene:
            last_gene_in_segment_feature = get_feature_by_gene_name(
                    last_gene, 'CDS', genome_record)
            last_gene_in_segment = get_feature_gene(
                    last_gene_in_segment_feature)
            hit_end_of_genome = _check_if_end_of_genome(last_gene_in_segment)
        else:
            last_gene_in_segment_feature = first_gene_in_segment_feature
            last_gene_in_segment = get_feature_gene(last_gene_in_segment_feature)
            seg_size = 0
            while seg_size < partition_size_lower_bound:
                assert last_gene_in_segment, "Unexpected CDS without gene name."
                # Step to next gene in segment and obtain the new segment size.
                # If next gene doesn't exist, just use what we have.
                hit_end_of_genome = _check_if_end_of_genome(
                        last_gene_in_segment)
                if hit_end_of_genome:
                    break
                last_gene_in_segment_feature = (
                        partition_helper.recoded_linked_cds_feature_map[
                                last_gene_in_segment])
                seg_size = (last_gene_in_segment_feature.location.end -
                        first_gene_in_segment_feature.location.start)
                last_gene_in_segment = get_feature_gene(
                        last_gene_in_segment_feature)

        # Now we look at the next few genes until we get a gene that starts a new
        # operon.
        if not should_heed_last_gene and not hit_end_of_genome:
            # TODO: Checking for hit_end_of_genome inside this block.
            first_gene_in_next_operon_feature = (
                    partition_helper.recoded_linked_cds_feature_map[
                            last_gene_in_segment])
            first_gene_in_next_operon = get_feature_gene(
                    first_gene_in_next_operon_feature)
            while not partition_helper.get_interval_before_gene(
                    first_gene_in_next_operon):
                # Shift the last segment, and next segment.
                last_gene_in_segment_feature = first_gene_in_next_operon_feature
                last_gene_in_segment = first_gene_in_next_operon
                first_gene_in_next_operon_feature = (
                    partition_helper.recoded_linked_cds_feature_map[
                            last_gene_in_segment])
                first_gene_in_next_operon = get_feature_gene(
                    first_gene_in_next_operon_feature)

        ### 2. Determine where to make the cut.

        # For now we choose the mid-point between the last gene in this
        # segment and the first gene in the next.
        # TODO: Account for promoter sites and other important motifs to
        #     make the cut.

        if hit_end_of_genome:
            recoded_seg_end = len(genome_record)
        else:
            recoded_seg_end = (
                    (first_gene_in_next_operon_feature.location.start +
                            last_gene_in_segment_feature.location.end) / 2)

        ### 3. Create annotation for the segment.
        seg_name = seg_annotation_prefix + str(start_seg_numbering + seg_index)
        seg_feature_location = FeatureLocation(recoded_seg_start, recoded_seg_end)
        add_feature_to_seq_record(genome_record, SeqFeature(
                seg_feature_location,
                type=InsertType.SYNTHESIS_SEG,
                id=seg_name,
                strand=1,
                qualifiers={'label': [seg_name]}
        ))

        # Add data to report.
        seg_size = recoded_seg_end - recoded_seg_start

        segment_report_list.append({
            'first_gene': first_gene_in_segment,
            'last_gene': last_gene_in_segment,
            'seg_start_bp': recoded_seg_start + 1,
            'seg_end_bp': recoded_seg_end,
            'seg_size': seg_size
        })

        if is_final_segment or hit_end_of_genome:
            break

        ### 4. Updates for next iteration of loop.
        first_gene_in_segment = first_gene_in_next_operon
        first_gene_in_segment_feature = first_gene_in_next_operon_feature
        mg1655_downstream_interval = (
                partition_helper.get_interval_before_gene(
                        first_gene_in_next_operon))
        mg1655_upstream_interval = mg1655_downstream_interval
        recoded_seg_start = recoded_seg_end

    # Write the output if requested.
    if annotated_genome_record_outfile:
        SeqIO.write(genome_record, annotated_genome_record_outfile, 'genbank')

    if output_segment_report:
        SEGMENT_REPORT_FIELDNAMES = [
            'first_gene',
            'last_gene',
            'seg_start_bp',
            'seg_end_bp',
            'seg_size'
        ]

        with open(output_segment_report, 'w') as report_fh:
            writer = csv.DictWriter(report_fh, SEGMENT_REPORT_FIELDNAMES)
            writer.writeheader()
            for segment in segment_report_list:
                writer.writerow(segment)


if __name__ == '__main__':
    import os
    from paths import DATA_DIR

    SEG5_2MB_DIR = os.path.join(DATA_DIR, 'completed_segments', 'seg5-seg44')

    SOURCE_RECORD_PATH = os.path.join(SEG5_2MB_DIR,
            '2013_08_27_13_16_44_mds42_refactored.gbk')

    genome_record = SeqIO.read(SOURCE_RECORD_PATH, 'genbank')

    ANNOTATED_OUTPUT = os.path.join(SEG5_2MB_DIR,
            '2013_08_28.seg5-seg44.genbank')

    OUTPUT_SEGMENT_REPORT = os.path.join(SEG5_2MB_DIR,
            '2013_08_28.seg5-seg44.segment_report.csv')

    find_segment_partitions(genome_record, num_segments=40, first_gene='yafD',
            start_seg_numbering=5, seg_annotation_prefix='seg',
            annotated_genome_record_outfile=ANNOTATED_OUTPUT,
            output_segment_report=OUTPUT_SEGMENT_REPORT)
