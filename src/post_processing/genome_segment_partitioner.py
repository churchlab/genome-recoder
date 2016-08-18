"""
Module that provides methods to determine the ideal Gen9 partitioning.

At a high-level, the strategy is to create some initial partition divider
objects that define fragments 3000 bases size, with 50 base overlaps for
assembly. We then iteratively check the underlying sequences as defined
by the partition dividers and wiggle their positions within some bounds
until we find a solution where overhang pieces have minimal secondary
structure and minimal homology with each other.
"""

import copy
import math
import sys

from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC

from biopython_util import add_feature_to_seq_record
from biopython_util import InsertType
from tool_wrappers.cd_hit_wrapper import CDHitWrapper
from tool_wrappers.hybrid_ss_min_wrapper import HybridSSMinWrapper

# TODO: Fix this for public release.
# # Append primer_design src to the path.
# sys.path.append('/home/glebk/Projects/churchlab/primer_design')
# from primer_design import find_primer


# The size of the overlap.
OVERHANG_SIZE = 50

# In order to avoid strong secondary structure, we want to get the
# ss value for the sequence to be between this and 0.
START_MIN_DELTA_G = -12

# The max allowed homology score.
START_MAX_PAIRWISE_HOMOLOGY = 0.6

# GC content bounds.
MIN_GC_CONTENT = 0.2
MAX_GC_CONTENT = 0.8

# Primer design
PRIMER_TARGET_TM = 65
MIN_PRIMER_SIZE = 17

# Subjective measure of score in our hand-rolled primer_design.score_primer()
# method. May need adjustment.
MIN_PRIMER_SCORE = 3

# We care about overlaps in terms of DNA.
EXTRA_UNAFOLD_OPTIONS = (
    '-n DNA -N 0.050 -M 0.0015 --mfold=50 --tmin=25 --tmax=25')


class MaxWiggleError(RuntimeError):
    """Exception thrown when the a divider has wiggled maximally."""
    pass


class NoValidPartitionError(RuntimeError):
    """Exception thrown when no valid partition solution is found."""
    pass



class PartitionedSequence(object):
    """Object representing a sequence partitioned with PartitionDivider
    objects.

    This object encapsulates most of the logic for finding the best
    partitioning solution.
    """
    # The underlying sequence.
    sequence = None

    # Stores a list of PartitionDivider objects.
    partition_divider_list = None

    def __init__(self, sequence, fragment_size, max_fragment_size, max_wiggle,
            minimize_fragments=True):
        """
        Args:
            sequence: Underlying sequence.
            fragment_size: Target fragment size.
            max_fragment_size: Maximum fragment size allowed.
        """
        # Set attributes.
        self.sequence = sequence
        self.fragment_size = fragment_size
        self.max_fragment_size = max_fragment_size
        self.max_wiggle = max_wiggle
        self.min_delta_g = START_MIN_DELTA_G
        self.max_homology = START_MAX_PAIRWISE_HOMOLOGY
        self.minimize_fragments = minimize_fragments

        # Testing min size.
        self.min_fragment_size = fragment_size - max_wiggle

        # Utility objects.
        self.ss_computer = HybridSSMinWrapper(extra_options=EXTRA_UNAFOLD_OPTIONS)
        cd_hit_wrapper_kwargs = {'threshold': self.max_homology}
        self.cd_hit_wrapper = CDHitWrapper(**cd_hit_wrapper_kwargs)

        # Methods to run after attributes are set.
        self.num_fragments = self.get_num_fragments_based_on_seq_length()


    def get_num_fragments_based_on_seq_length(self):
        """Figure out how many fragments the given sequence will be
        partitioned into.

        NOTE: This currently assumes that this wont change during the algo
            (i.e. we always have 1 short piece and the rest are fragment_size).
        """
        normal_fragment_size_segs = int(math.ceil(float(len(self.sequence)) /
                (self.fragment_size - OVERHANG_SIZE)))

        if not self.minimize_fragments:
            return normal_fragment_size_segs

        # See if we can do 1 less fragment by getting close to the max
        # fragment size. Arbitrarily OVERHANG_SIZE less than max.
        squeeze_fragment_size = self.max_fragment_size - OVERHANG_SIZE
        squeeze_num_segs = int(math.ceil(float(len(self.sequence)) /
                (squeeze_fragment_size - OVERHANG_SIZE)))

        if squeeze_num_segs < normal_fragment_size_segs:
            self.fragment_size = squeeze_fragment_size
            self.max_wiggle = self.max_fragment_size - self.fragment_size
            return squeeze_num_segs

        # Otherwise just return the normal number.
        return normal_fragment_size_segs


    def get_fragment_list(self):
        """Method that does the actual work of finding the fragments.
        """
        # Iterate over all the possible places to put the small leftover fragment.
        # This gives us a bit more flexibility in finding a possible solution,
        # without needing to go through the effort of setting up a proper
        # graph search.
        for leftover_fragment_index in range(self.num_fragments):
            print 'Iterating on leftover_fragment_index %d of %d' % (
                    leftover_fragment_index + 1, self.num_fragments)
            # print '...leftover_fragment_index %d' % leftover_fragment_index
            try:
                # Generate the PartitionDividers with the leftover fragment
                # having the currently iterated index.
                self._generate_partition_divider_list(leftover_fragment_index)

                # Search for the optimal partitioning.
                self._find_optimal_partitioning()

                # Now that we've found an optimal partitioning, wiggle a bit
                # more to optimize for primer design.
                # NOTE: Be careful to avoid going beyond max fragment size at
                # this step.
                self._optimize_partitions_for_primer_design()

            except MaxWiggleError:
                continue


            # Get the resulting fragments.
            fragment_seq_list = self._get_current_fragment_list()

            # Check the fragment list.
            self._assert_current_fragment_reassembly(fragment_seq_list)

            return fragment_seq_list


        # TODO: Rather than just failing, bump the threshold depending on the
        # threshold that caused the most errors in this run.

        raise NoValidPartitionError()


    def _generate_partition_divider_list(self, leftover_fragment_index):
        """Reset the list of PartitionDivider objects, with the small
        left-over fragment placed at the index specified.

        This lets us explore alternate possibilities for partitioning. The
        resulting partitioning may not be optimal. To find the optimal
        paritioning, follow up a call to this method with a call to
        _find_optimal_partitioning().

        Args:
            leftover_fragment_index: Index at which to insert the leftover
                fragment, that is the typically < 3kb piece (the current
                GENE_BYTE) size as of writing this code).
        """
        self.partition_divider_list = []

        total_fragments = self.get_num_fragments_based_on_seq_length()

        # We want to figure out the size of the leftover bit that won't be
        # enough for a full fragment (e.g. 3000 bases for Gen9).
        # Each fragment includes the overhang from the next
        # one, so to avoid double-counting, we consider each fragment to be
        # fragment_size - OVERHANG_SIZE (e.g. 2950 bases for Gen9).
        leftover_fragment_size = max(
                self.min_fragment_size,
                len(self.sequence) - (total_fragments - 1) *
                        (self.fragment_size - OVERHANG_SIZE))
        assert leftover_fragment_size < self.max_fragment_size, (
                "Actual leftover fragment size: %d" % leftover_fragment_size)

        # Adjust the fragment size given the new leftover fragment size.
        self.fragment_size = (
                (len(self.sequence) - leftover_fragment_size) /
                        (total_fragments - 1) + OVERHANG_SIZE)

        # Figure out the positions for the dividers, with attention to where
        # to put the short fragment.
        current_fragment_index = 0
        next_start_position = 0
        while current_fragment_index < total_fragments - 1:
            if current_fragment_index == leftover_fragment_index:
                next_start_position += leftover_fragment_size - OVERHANG_SIZE
            else:
                next_start_position += self.fragment_size - OVERHANG_SIZE
            self.partition_divider_list.append(
                    PartitionDivider(next_start_position, self.max_wiggle))
            current_fragment_index += 1


    def get_partition_divider_list(self):
        """Returns the list of PartitionDivider objects.
        """
        return self.partition_divider_list


    def get_overhang_sequence_list(self):
        """Returns the list of overhang sequences as defined by the current
        PartitionDivider settings.
        """
        return map(self.get_seq_for_partition_divider,
            self.partition_divider_list
        )


    def get_seq_for_partition_divider(self, partition_divider):
        """Returns the portion of sequence corresponding to the current
        PartitionDivider state.
        """
        return self.sequence[partition_divider.get_start():
                partition_divider.get_end()]


    def _get_current_fragment_list(self):
        """Returns the list of Gen9 fragments as defined by the current
        state of the PartitionDividers.
        """
        fragment_seq_list = []

        # Do the first one manually.
        first_divider = self.partition_divider_list[0]
        first_seq = self.sequence[:first_divider.get_end()]
        fragment_seq_list.append({
                'sequence': first_seq,
                'start': 0,
                'end': len(first_seq)
        })

        # The rest we do as defined by the current divider as the
        # starting point and the next one as the end point.
        for divider_index in range(len(self.partition_divider_list) - 1):
            start_divider = self.partition_divider_list[divider_index]
            end_divider = self.partition_divider_list[divider_index + 1]
            fragment_seq_list.append({
                'sequence': self.sequence[
                        start_divider.get_start():end_divider.get_end()],
                'start': start_divider.get_start(),
                'end': end_divider.get_end()
            })

        # And do the last one manually.
        last_divider = self.partition_divider_list[-1]
        last_seq = self.sequence[last_divider.get_start():]
        fragment_seq_list.append({
                'sequence': last_seq,
                'start': last_divider.get_start(),
                'end': len(self.sequence)
        })

        assert len(self.partition_divider_list) + 1 == len(fragment_seq_list)

        return fragment_seq_list


    def _assert_current_fragment_reassembly(self, fragment_seq_list=None):
        """Make sure the original genome sequence by reassembling the fragments
        with overlaps of OVERHANG_SIZE.
        """
        if not fragment_seq_list:
            fragment_seq_list = self._get_current_fragment_list()

        # Check fragment lengths are acceptable.
        fragment_seq_lengths = [len(fragment['sequence']) for fragment in
                fragment_seq_list]
        for i, frag_len in enumerate(fragment_seq_lengths):
            assert frag_len < self.max_fragment_size, (
                    "Segment index %d failed with length %d" % (i, frag_len))

        # Check that reassembling the fragments, accounting for overhangs,
        # is correct.
        whole_assembly = self._assemble_sequence_from_fragments(
                fragment_seq_list)

        assert self.sequence == whole_assembly, (
                "Re-assembled sequence doesn't match original.")


    def _assemble_sequence_from_fragments(self, fragment_seq_list):
        """Helper method for reassembling fragments into a sequence.

        Mostly used for assertions/testing.
        """
        # Check that reassembling the fragments, accounting for overhangs,
        # is correct.
        whole_assembly = ''

        # Treat the first one differently.
        whole_assembly += fragment_seq_list[0]['sequence']

        for idx in range(1, len(fragment_seq_list)):
            whole_assembly += fragment_seq_list[idx]['sequence'][OVERHANG_SIZE:]

        return whole_assembly


    def _find_optimal_partitioning(self):
        """Searches for the optimal partioning.
        """
        # Hack to allow iterating on troubled index first.
        hack_granular_wiggle_trouble_idx = 0

        # Repeat until no more change are made.
        at_least_one_wiggle = True
        while at_least_one_wiggle:
            at_least_one_wiggle = False
            at_least_one_wiggle = (at_least_one_wiggle or
                    self._resolve_independent_constraint_violations())
            at_least_one_wiggle = (at_least_one_wiggle or
                    self._resolve_homology_issues())
            # # This next check is expensive so don't do it unless
            # # above checks did not require wiggle on this iteration.
            if at_least_one_wiggle:
                continue
            at_least_one_wiggle, hack_granular_wiggle_trouble_idx = (
                    self._granular_wiggle_for_good_primers(
                            do_first=hack_granular_wiggle_trouble_idx))


    def _resolve_independent_constraint_violations(self):
        """Loops through the PartitionDivider objects, wiggling them
        in order to resolve any constraints violations.

        Returns:
            Boolean indicating whether there was at least one wiggle.
        """
        at_least_one_wiggle = False
        for partition_divider in self.partition_divider_list:
            # Loop to make sure fixing one doesn't break the other.
            wiggled = True
            while wiggled:
                wiggled = False
                wiggled = wiggled or self._resolve_secondary_structure(
                        partition_divider)
                wiggled = wiggled or self._resolve_gc_content(
                        partition_divider)
                if wiggled:
                    at_least_one_wiggle = True
        return at_least_one_wiggle


    def _resolve_secondary_structure(self, partition_divider):
        """Checks that the overhang sequence represented by a PartitionDivider
        has minimal secondary structure as approximated by delta_g, wiggling
        as necessary to get it to fit.

        Args:
            partition_divider: A PartitionDivider object that corresponds
                to an overhang sequence.

        Returns:
            Boolean indicating whether a wiggle took place.
        """
        wiggled = False
        overhang_seq = self.get_seq_for_partition_divider(partition_divider)
        delta_g = self.ss_computer.compute_delta_g(overhang_seq)
        while delta_g < self.min_delta_g:
            partition_divider.wiggle()
            wiggled = True
            overhang_seq = self.get_seq_for_partition_divider(partition_divider)
            delta_g = self.ss_computer.compute_delta_g(overhang_seq)
        return wiggled


    def _resolve_gc_content(self, partition_divider):
        """Checks that the overhang sequence represented by a PartitionDivider
        has a GC content within the desired range, wiggling as necessary.

        Args:
            partition_divider: A PartitionDivider object that corresponds
                to an overhang sequence.

        Returns:
            Boolean indicating whether a wiggle took place.

        """
        wiggled = False
        overhang_seq = self.get_seq_for_partition_divider(partition_divider)
        gc_content = GC(overhang_seq) / 100
        while gc_content < MIN_GC_CONTENT or gc_content > MAX_GC_CONTENT:
            partition_divider.wiggle()
            wiggled = True
            overhang_seq = self.get_seq_for_partition_divider(partition_divider)
            gc_content = GC(overhang_seq) / 100
        return wiggled


    def _resolve_homology_issues(self):
        """Wiggles until homology issues are resolved.

        Returns:
            Boolean indicating whether there was at least one wiggle.
        """
        at_least_one_wiggle = False
        has_homology_conflicts = True
        while has_homology_conflicts:
            has_homology_conflicts = False
            overhang_seq_list = self.get_overhang_sequence_list()
            conflict_set_list = self.cd_hit_wrapper.detect_high_homology_pairs(
                    overhang_seq_list)
            if len(conflict_set_list) > 0:
                at_least_one_wiggle = True
                has_homology_conflicts = True
                overhang_seq_list = self.get_overhang_sequence_list()
                for conflict_set in conflict_set_list:
                    for idx in conflict_set:
                        self.partition_divider_list[idx].wiggle()
        return at_least_one_wiggle


    def _granular_wiggle_for_good_primers(self, do_first=0):
        """Small wiggles to make sure we can design good primers for the
        fragments.

        For each partition, try to find good primers. If not found,
        do a big wiggle and report at_least_one_wiggle.

        To find the best primers, we wiggle each partition divider over a small
        range 2 * EPSILON, trying to find the best combo of fwd and reverse
        primers determined by that partition divider.

        Returns:
            A tuple pair, with the following ordered elements:
                * boolean of whether we did a big wiggle
                * int representing trouble index, 0 if none

        """
        EPSILON = 10

        def _find_primers_for_partition_divider(pd_idx, partition_divider):
            """Helper function that does the work in the loop that follows.
            """
            initial_start = partition_divider.get_start()
            start_range = range(initial_start - EPSILON,
                    initial_start + EPSILON)
            best_score = None
            for start in start_range:
                partition_divider.set_start(start)
                overhang_seq = str(self.get_seq_for_partition_divider(
                        partition_divider))

                f_pr_seq, f_pr_tm, f_pr_ss, f_pr_score = find_primer(
                            overhang_seq,
                            min_size=MIN_PRIMER_SIZE,
                            target_tm=PRIMER_TARGET_TM)

                r_pr_seq, r_pr_tm, r_pr_ss, r_pr_score = find_primer(
                            reverse_complement(overhang_seq),
                            min_size=MIN_PRIMER_SIZE,
                            target_tm=PRIMER_TARGET_TM)

                # We want a combined score function that gives a bonus
                # when both scores are medium high. In other words,
                # we want to avoid the case where one of the scores is
                # really good but the other one is not necessarily that good,
                # but is masked by the goodness of the other.
                combined_score = min(f_pr_score, r_pr_score)
                if not best_score or combined_score > best_score:
                    best_score = combined_score
                    best_start = start

            # Check if the best score is good enough.
            if best_score >= MIN_PRIMER_SCORE:
                partition_divider.set_start(best_start)
                return False, 0
            else:
                partition_divider.set_start(initial_start)
                partition_divider.wiggle()
                return True, pd_idx

        # Count of total divider for print statements.
        total_partition_dividers = len(self.partition_divider_list)

        # Hack to allow repeating the trouble one first.
        if do_first > 0:
            print 'optimizing %d of %d partitions' % (
                    do_first + 1, total_partition_dividers)
            partition_divider = self.partition_divider_list[do_first]
            at_least_one_wiggle, trouble_index = (
                    _find_primers_for_partition_divider(
                            do_first, partition_divider))
            if at_least_one_wiggle:
                return at_least_one_wiggle, trouble_index

        # Must pass for all in order to succeed.
        for idx, partition_divider in enumerate(self.partition_divider_list):
            print 'optimizing %d of %d partitions' % (
                    idx + 1, total_partition_dividers)
            at_least_one_wiggle, trouble_index = (
                    _find_primers_for_partition_divider(
                            idx, partition_divider))
            if at_least_one_wiggle:
                return at_least_one_wiggle, trouble_index

        # All succeeded.
        return False, 0


    def _optimize_partitions_for_primer_design(self):
        """Move the partitions around a bit to allow for better primers.

        raises:
            MaxWiggleError if valid primers not found.
        """
        EPSILON = 10
        # Try to make the best primers for the partition.
        # Wiggle each partition divider a small amount to find the
        # best combo of primers (one fwd one reverse), e.g.:
        #     |>>>>>>>.....<<<<<<|
        total_partition_dividers = len(self.partition_divider_list)
        for idx, partition_divider in enumerate(self.partition_divider_list):
            print 'optimizing %d of %d partitions' % (
                    idx + 1, total_partition_dividers)
            initial_start = partition_divider.get_start()
            start_range = range(initial_start - EPSILON,
                    initial_start + EPSILON)
            best_score = None
            for start in start_range:
                partition_divider.set_start(start)
                overhang_seq = self.get_seq_for_partition_divider(
                        partition_divider)

                f_pr_seq, f_pr_tm, f_pr_ss, f_pr_score = find_primer(
                            str(overhang_seq),
                            min_size=MIN_PRIMER_SIZE,
                            target_tm=PRIMER_TARGET_TM)

                r_pr_seq, r_pr_tm, r_pr_ss, r_pr_score = find_primer(
                            reverse_complement(str(overhang_seq)),
                            min_size=MIN_PRIMER_SIZE,
                            target_tm=PRIMER_TARGET_TM)

                # We want a combined score function that gives a bonus
                # when both scores are medium high. In other words,
                # we want to avoid the case where one of the scores is
                # really good but the other one is not necessarily that good,
                # but is masked by the goodness of the other.
                # combined_score = (f_pr_score + r_pr_score) ** 2
                combined_score = min(f_pr_score, r_pr_score)
                if not best_score or combined_score > best_score:
                    best_score = combined_score
                    best_start = start
            # Now that we've tried all start positions, check whether the best
            # one is good enough.
            if best_score < MIN_PRIMER_SCORE:
                print best_score
                raise MaxWiggleError()
            partition_divider.set_start(best_start)


class EvenlyPartitionedSequence(PartitionedSequence):
    """Represents a partitoined sequence that breaks the given segment as
    evenly as possible into the given number of fragments.

    This differs from PartitionedSequence which takes a set size, and then
    catches the leftover sequence into a remaining fragment.

    NOTE: I'm not sure whether sub-classing PartitionedSequence is the best way
    to implement this at the moment.
    """

    def __init__(self, sequence, num_fragments, overhang_size=50,
            max_wiggle=100, min_delta_g=-10):
        # Passed-in attributes.
        self.sequence = sequence
        self.num_fragments = num_fragments
        self.overhang_size = overhang_size
        self.max_wiggle = max_wiggle
        self.partition_divider_list = None
        self.min_delta_g = min_delta_g
        self.max_homology = START_MAX_PAIRWISE_HOMOLOGY

        # Calculated attributes.
        self.base_fragment_size = len(self.sequence) / self.num_fragments

        # Utility objects.
        self.ss_computer = HybridSSMinWrapper(extra_options=EXTRA_UNAFOLD_OPTIONS)
        cd_hit_wrapper_kwargs = {'threshold': self.max_homology}
        self.cd_hit_wrapper = CDHitWrapper(**cd_hit_wrapper_kwargs)

    def get_fragment_list(self):
        try:
            self._generate_partition_divider_list()
            self._find_optimal_partitioning()
            fragment_seq_list = self._get_current_fragment_list()
            return fragment_seq_list
        except MaxWiggleError:
            raise NoValidPartitionError()

    def _generate_partition_divider_list(self):
        """Partitions approximately evenly.
        """
        self.partition_divider_list = []
        next_partition = 0
        for i in range(self.num_fragments - 1):
            next_partition += self.base_fragment_size
            partition_divider_start = next_partition - self.overhang_size / 2
            self.partition_divider_list.append(
                    PartitionDivider(partition_divider_start, self.max_wiggle,
                    self.overhang_size))


class PartitionDivider(object):
    """Object with a root position and the ability to be slid with some
    flexibility to alleviate secondary structure or homology issues.

    This object effectively embodies the 50-base overhang shared by
    two adjacent fragments of the total 48-kbase sequence.
    """
    # The initial position of the divider.
    base_position = None

    # State variable that indicates the next wiggle to try when the overhang
    # sequence associated with this divider fails either a secondary structure
    # check or a homology check. Without loss of generality, the algorithm for
    # updating this will be to just negate it if it's negative and increment
    # by next_wiggle_inc and negate when positive, until the max is reached.
    next_wiggle_delta = -10

    # The amount to increment the wiggle.
    wiggle_inc = 10

    def __init__(self, base_position, max_wiggle, size=OVERHANG_SIZE):
        """Construct the object with its current_position set to base_position.
        """
        self.base_position = base_position
        self.max_wiggle = max_wiggle
        self.size = size

        self.current_start = self.base_position


    def __repr__(self):
        return 'PartitionDivider(base: %s, current %s)' % (
                self.base_position, self.current_start)


    def get_start(self):
        """Returns the current start position of the divider.
        """
        return self.current_start


    def set_start(self, start):
        self.current_start = start


    def get_end(self):
        """Returns the current end position of the divider.
        """
        return self.current_start + self.size


    def wiggle(self):
        """Wiggle the position of the divider to try a different partition.
        Used when resolving a secondary structure or homology issue.
        """
        if abs(self.next_wiggle_delta) > self.max_wiggle:
            raise MaxWiggleError()

        self.current_start = self.base_position + self.next_wiggle_delta
        self._update_wiggle()


    def _update_wiggle(self):
        """Updates the next wiggle delta.

        This is done by negating the current value, and also incrementing
        the absolute value on positive cycles (without loss of generality).
        """
        if self.next_wiggle_delta > 0:
            self.next_wiggle_delta += self.wiggle_inc
        self.next_wiggle_delta = -1 * self.next_wiggle_delta



def partition_segment(
        genome_record,
        start_position,
        end_position,
        fragment_size,
        max_fragment_size,
        max_wiggle,
        annotated_genome_record_outfile,
        fragments_outfile,
        fragment_annotation_id_prefix,
        validation_start_seq=None,
        validation_end_seq=None,
        mutate_genome_record=False):
    """Given a segment of size 48kb, breaks it up into chunks of size 3kb
    and prepares the order to be sent to Gen9. The resulting fragments to be
    synthesized is written to the output file.

    Visually the output would look something like this:

        >>>>>>>>>>>>>>>>>>>>>>>>>
        hhhh|               |tttt

                             >>>>>>>>>>>>>>>>>>>>>>>
                             hhhh|             |tttt

                                                >>>>>>>>>>>>>>>>>>>>>>>
                                                hhhh|             |tttt

    Args:
        genome_record: Source SeqRecord for the segment to be created.
        start_position: Start position relative to genome_record for the entire
            segment.
        end_position: End position relative to genome_record for the entire
            segment.
        annotated_genome_record_outfile: File to which the annotated record
            will be written to. The annotations will be for each separate
            gene piece to be synthesized by Gen9. Pass None to skip.
        fragments_outfile: Outfile for the gene pieces to be sent to Gen9.
            Pass None to skip writing outfile.
        fragment_annotation_id_prefix: Prefix for annotations added.
            (e.g 'seg2')
        validation_start_seq: Optional sequence at the start of the segment
            to sanity check the positions.
        validation_end_seq: Optional sequence at the end of the segment
            to sanity check the positions.
        mutate_genome_record: Modify the genome record in-place rather than
            making a copy.
    """
    # HACK: seg4 - seg48 order, seg25 was hitting 4200 bound. Just do it the
    # dumb way.
    minimize_fragments = True
    if fragment_annotation_id_prefix == 'seg25':
        minimize_fragments = False

    orig_seq = str(genome_record.seq)

    # Maybe perform validation.
    if validation_start_seq:
        actual_start_seq = orig_seq[
                start_position:start_position + len(validation_start_seq)]
        assert validation_start_seq == actual_start_seq, (
                "Actual: %s" % actual_start_seq)

    if validation_end_seq:
        actual_end_seq = orig_seq[
                end_position - len(validation_end_seq):end_position]
        assert validation_end_seq == actual_end_seq, (
                "Actual: %s" % actual_end_seq)


    ### Find a valid partitioning.

    # The strategy is to determine an initial set of PartitionDivider objects,
    # each of which has a starting position, and then can be slide back and
    # forth some amount if necessary to alleviate secondary structure or
    # homology issues. The PartitionedSequence object encapsulates the logic
    # to do this.
    forward_seq = str(genome_record.seq[start_position:end_position])
    partitioned_seq_obj = PartitionedSequence(
            forward_seq, fragment_size, max_fragment_size, max_wiggle,
            minimize_fragments=minimize_fragments)
    fragment_list = partitioned_seq_obj.get_fragment_list()


    ### Write the results.

    # Turn the objects into SeqRecords to give them ids.
    fragments_as_seq_records = []
    fragment_features = []
    for idx, fragment_obj in enumerate(fragment_list):
        fragment_str = fragment_obj['sequence']
        fragment_seq = Seq(fragment_str, generic_dna)

        # Generate a somewhat meaningful id.
        fragment_idx_sring = ('00' + str(idx))[-3:]
        fragment_id = fragment_annotation_id_prefix + '_' + fragment_idx_sring

        # And some brief info in the description
        fragment_description = "Size: %d" % len(fragment_str)

        # Append to the list that we'll write out as .fasta.
        fragments_as_seq_records.append(
                SeqRecord(fragment_seq, id=fragment_id,
                        description=fragment_description))

        # Create a SeqFeature that we'll add to the annotation out file.
        fragment_feature = SeqFeature(
                location=FeatureLocation(
                        fragment_obj['start'] + start_position,
                        fragment_obj['end'] + start_position),
                type=InsertType.PARTITION_GEN9_SEG,
                strand=1,
                id=fragment_id)
        fragment_feature.qualifiers['label'] = [fragment_id]
        fragment_features.append(fragment_feature)

    if fragments_outfile is not None:
        with open(fragments_outfile, 'w') as fh:
            SeqIO.write(fragments_as_seq_records, fh, 'fasta')

    # Add feature annotations to a copy of genome record and write to the
    # outfile.
    if mutate_genome_record:
        genome_record_copy = genome_record
    else:
        genome_record_copy = copy.deepcopy(genome_record)

    for fragment_feature in fragment_features:
        add_feature_to_seq_record(genome_record_copy, fragment_feature)

    if annotated_genome_record_outfile is not None:
        with open(annotated_genome_record_outfile, 'w') as fh:
            SeqIO.write(genome_record_copy, fh, 'genbank')

    print '...Done.'


if __name__ == '__main__':
    from biopython_util import get_genome_record

    GENOME_RECORD = get_genome_record(
            '../data/completed_segments/seg2/seg2_final_with_flanking.gbk')

    START_POSITION = 100863
    END_POSITION = 149941

    FRAGMENT_SIZE = 3000
    MAX_FRAGMENT_SIZE = 3100

    VALIDATION_START_SEQ = 'CAGCCTTGTTTCGCCA'
    VALIDATION_END_SEQ = 'AGAATCATGGCCT'

    ANNOTATED_GENOME_RECORD_OUTFILE = 'output/test_gen9_out.gbk'
    FRAGMENTS_OUTFILE = 'output/seg2_gen9_fragments.fa'
    FRAGMENT_ANNOTATION_ID_PREFIX = 'seg2_'

    partition_segment(
            GENOME_RECORD,
            START_POSITION,
            END_POSITION,
            FRAGMENT_SIZE,
            MAX_FRAGMENT_SIZE,
            ANNOTATED_GENOME_RECORD_OUTFILE,
            FRAGMENTS_OUTFILE,
            FRAGMENT_ANNOTATION_ID_PREFIX,
            VALIDATION_START_SEQ,
            VALIDATION_END_SEQ
    )
