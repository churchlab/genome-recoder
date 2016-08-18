"""
Methods for fixing the GC content of a Genome based on constraints imposed by,
say, synthesis requirements.
"""

import copy
import csv
import random

from Bio.SeqUtils import GC

from biopython_util import calc_interval_list_to_features_overlapped
from biopython_util import get_region_codon_indeces_in_feature
from biopython_util import update_seq_record_feature
from codon_replacer_util import replace_codons_in_single_feature
from genome_diff import find_mismatches_between_same_size_genomes


GC_VALUES = set(['G', 'C', 'g', 'c'])


def get_GC_optimized(seq, start, end, GC_one_left=None):
    """Calculates GC content for the sequence in place. Actually,
    in-place means rather than passing the sliced sequence, you pass
    the start and end positions in the given seq.

    Args:
        seq: The sequence.
        start: The start index into the sequence.
        end: The end index into the sequence.
        GC_one_left: If provided, implies:
            * Last interval calculated was one to the left.
            * Interval size is the same.

    Returns:
        Float between 0 and 1.
    """
    interval_size = end - start
    if GC_one_left is not None:
        GC_count = int(GC_one_left * interval_size)
        if seq[start - 1] in GC_VALUES:
            GC_count -= 1
        if seq[end] in GC_VALUES:
            GC_count += 1
    else:
        GC_count = 0
        for pos in xrange(start, end):
            if seq[pos] in GC_VALUES:
                GC_count += 1
    return GC_count / interval_size


class GCContentConstraints(object):
    """Object that defines the constraints for GC content.

    Clients should construct one of these and pass it to the fixer.
    This is better than having a plethora of params.
    """
    # Defaults based on Gen9 Faq:
    DEFAULT_LOCAL_WINDOW_SIZE = 100
    DEFAULT_LOCAL_LOWER_BOUND = 0.30
    DEFAULT_LOCAL_UPPER_BOUND = 0.75

    def __init__(self):
        self.local_window_size = self.DEFAULT_LOCAL_WINDOW_SIZE
        self.local_window_lower_bound = self.DEFAULT_LOCAL_LOWER_BOUND
        self.local_window_upper_bound = self.DEFAULT_LOCAL_UPPER_BOUND

    # TODO: Add setters and getters as needed.


def find_gc_content_extremes(
        genome_record,
        gc_content_constraint_obj=GCContentConstraints(),
        start_bound=None,
        end_bound=None,
        debug=False):
    """Finds runs of extreme GC content.

    Args:
        genome_record: The SeqRecord object with the sequence.
        gc_content_constraint_obj: A GCContentConstraints object that
            allows the client to configure the fixes.
        start_bound: Optionally bound fixes to start at this position.
        end_bound: Optionally bound fixes to end at this position.

    Returns:
        List of objects with keys:
            * interval: Pythonic (start, end) of the interval.
            * avg_gc: Average GC content over this interval.
    """
    extreme_gc_intervals = []

    effective_start_bound = start_bound if start_bound else 0
    effective_end_bound = end_bound if end_bound else len(genome_record.seq)

    running_interval = None
    running_gc_total = 0

    # Slide the window looking for violations of GC content restrictions.
    window_center_range = xrange(
            effective_start_bound +
                    gc_content_constraint_obj.local_window_size / 2,
            effective_end_bound -
                    gc_content_constraint_obj.local_window_size / 2)
    # Necessary initialization for our get_GC_optimized() method.
    gc_content = None
    for window_center_pos in window_center_range:
        window_start_pos = (window_center_pos -
                gc_content_constraint_obj.local_window_size / 2)
        window_end_pos = (window_start_pos +
                gc_content_constraint_obj.local_window_size)
        gc_content = GC(genome_record.seq, window_start_pos,
                window_end_pos, gc_content)
        if (gc_content_constraint_obj.local_window_lower_bound <= gc_content <=
                gc_content_constraint_obj.local_window_upper_bound):
            # End of extreme interval. Record the current interval and reset.
            if running_interval:
                interval_size = running_interval[1] - running_interval[0] + 1
                avg_gc = running_gc_total / interval_size
                extreme_gc_intervals.append({
                    'interval': running_interval,
                    'avg_gc': avg_gc
                })

                # Reset.
                running_interval = None
                running_gc_total = 0
        else:
            # Create or update the running interval.
            if not running_interval:
                running_interval =  (window_center_pos, window_center_pos)
            else:
                running_interval =  (running_interval[0], window_center_pos)
            running_gc_total += gc_content

    return extreme_gc_intervals


def fix_gc_content(
        refactor_context,
        gc_content_constraint_obj,
        start_bound=None,
        end_bound=None,
        debug=False,
        report_file=None):
    """Fixes the GC content according to desired constraints.

    Strategy:
        Slide a window across the genome and bump any regions that fall
        outside of the constraint. Now, there is some subtlety here in that
        we have a good idea of how to fix coding regions (i.e. synonymous
        codon swaps), but want to avoid messing with stuff outside of
        coding regions. And so when we identify a bad window, for now,
        we limit fixes to any coding portions of that window only.

    TODOs:
        * We are only dealing with local window for now. Figure out how we want
        to deal with global window.

    Args:
        refactor_context: The RefactorContext.
        gc_content_constraint_obj: A GCContentConstraints object that
            allows the client to configure the fixes.
        start_bound: Optionally bound fixes to start at this position.
        end_bound: Optionally bound fixes to end at this position.
        debug: Debug flag. Prints helpful output. For now, runs analysis only,
            and doesn't actually make changes.

    Returns:
        A copy of the genome_record contained within refactor_context
        with the GC content made to satisfy constraints.
    """
    print 'Fixing GC content...'
    updated_genome_record = copy.deepcopy(refactor_context.get_genome_record())

    # Figure out effective bounds.
    effective_start_bound = start_bound if start_bound else 0
    effective_end_bound = end_bound if end_bound else len(updated_genome_record)

    # Features that we can do synonymous swaps in
    swappable_features = [feature for feature in updated_genome_record.features
            if feature.type == 'CDS']

    # Slide the window looking for violations of GC content restrictions.
    window_center_range = range(
            effective_start_bound +
                    gc_content_constraint_obj.local_window_size / 2,
            effective_end_bound -
                    gc_content_constraint_obj.local_window_size / 2)
    report_intervals = []
    if debug:
        running_interval = None
        running_gc_total = 0
    for window_center_pos in window_center_range:
        window_start_pos = (window_center_pos -
                gc_content_constraint_obj.local_window_size / 2)
        window_end_pos = (window_start_pos +
                gc_content_constraint_obj.local_window_size)
        window_seq = updated_genome_record.seq[window_start_pos:window_end_pos]
        gc_content = GC(window_seq) / 100
        if (gc_content_constraint_obj.local_window_lower_bound <= gc_content <=
                gc_content_constraint_obj.local_window_upper_bound):
            # GC is all good.
            if debug:
                # Close the running interval and print it out.
                if running_interval:
                    interval_size = running_interval[1] - running_interval[0] + 1
                    avg_gc = running_gc_total / interval_size
                    report_intervals.append({
                        'interval': str(running_interval),
                        'interval_size': interval_size,
                        'avg_gc': avg_gc,
                    })
                    print ('%s, size: %d, average_gc: %f' % (
                            str(running_interval),
                            interval_size,
                            avg_gc))
                    running_interval = None
                    running_gc_total = 0
            continue

        if debug:
            if not running_interval:
                running_interval =  (window_center_pos, window_center_pos)
            else:
                running_interval =  (running_interval[0], window_center_pos)
            running_gc_total += gc_content
            continue

        # As a first stab, only attempt fixes in the simplest of cases.
        # That is, only do synonymous codon swaps within parts of features
        # that are not overlapping.

        # First identify all features overlaped by the interval.
        interval = (window_start_pos, window_end_pos)
        overlapped_features = calc_interval_list_to_features_overlapped(
                [interval], swappable_features)[0]
        if len(overlapped_features) != 1:
            # TODO: Eventually handle more complex cases.
            continue

        # Otherwise attempt to fix.
        feature = overlapped_features[0]
        feature_seq = str(feature.extract(updated_genome_record.seq))

        # Figure out the specific codons that need to be changed.
        affected_codon_indeces = get_region_codon_indeces_in_feature(
                feature, interval)
        avoid_codons_in_positions = {}
        for codon_index in affected_codon_indeces:
            codon = feature_seq[codon_index * 3 : codon_index * 3 + 3]
            if GC(codon) < 1.0:
                avoid_codons_in_positions[codon_index] = codon

        # Perform replace.
        first_codon_to_modify = affected_codon_indeces[0]
        last_codon_to_modify = affected_codon_indeces[-1]
        assert first_codon_to_modify <= last_codon_to_modify
        result = replace_codons_in_single_feature(
                refactor_context,
                feature.id,
                explicit_genome_record=updated_genome_record,
                start_codon_index=first_codon_to_modify,
                last_codon_index=last_codon_to_modify,
                avoid_codons_in_positions=avoid_codons_in_positions)
        if not result['is_success']:
            # TODO: Do something better for debugging here, although
            # we don't necessarily need each replace to succeed.
            continue

        update_seq_record_feature(
                updated_genome_record,
                feature.id,
                result
        )

    print '...Done.'

    if report_file:
        print 'Writing report.'
        REPORT_FIELDNAMES = [
            'interval',
            'interval_size',
            'avg_gc',
        ]
        with open(report_file, 'w') as report_fh:
            writer = csv.DictWriter(report_fh, REPORT_FIELDNAMES)
            writer.writeheader()
            for interval in report_intervals:
                writer.writerow(interval)

    return updated_genome_record


def automated_intergenic_gc_fixer(genome_record, interval_list,
        gc_content_constraint_obj=GCContentConstraints(),
        changes_ok_in_feature_types=[], increase_GC=True):
    """Fixes GC content in the given interval.

    Avoids making changes within any feature annotations and 20 bp upstream
    of CDS, unless feature type is specified in changes_ok_in_feature_types.

    Args:
        genome_record: Mutable SeqRecord.
        interval_list: List of tuples in which we want to fix GC. Note that
            bases just outside this interval may be changed, as they
            contribute to the local GC measure.
        changes_ok_in_feature_types: List of feature types that we want to
            allow changes inside of.
        gc_content_constraint_obj: Provide limits on GC content.
        increase_GC: If True, increase GC. If False, decrease GC.

    Limitations:
        * Only changes bases that are not inside of any feature annotation
          or before CDS gene.
        * Only increases or decreases GC in all given intervals.
    """
    # Strategy: Sample positions in the interval and change them to
    # corresponding purine or pyrimidine. Recalculate interval GC until
    # above threshold.

    assert increase_GC, "Implementation only supports increasing GC right now."

    # Record the original length for a final assertion.
    len_genome_before_fix = len(genome_record)

    # First calculate black-listed positions that cannot be changed. These are
    # positions that are inside of feature annotations, or just upstream of CDS.
    pos_blacklist_set = _calculate_feature_annotation_shadow(genome_record,
            changes_ok_in_feature_types=changes_ok_in_feature_types)

    # We extend each target interval by half of the GC measurement window size
    # in each direction, since the we measure the GC centered at each position
    # in the interval.
    half_window = gc_content_constraint_obj.local_window_size / 2

    # R -> R, Y -> Y
    GC_TRANSITION_TABLE = {
        'A': 'G',
        'T': 'C',
        'G': 'A',
        'C': 'T',
    }

    GC_bases = 'GC'

    for interval in interval_list:
        # Extend the interval to include the epsilon on each side. We change
        # positions in this extended_interval, although we still only target
        # raising the GC of windows centered about positions in the unextended
        # interval.
        extended_interval = (interval[0] - half_window,
                interval[1] + half_window)
        interval_size = interval[1] - interval[0]
        extended_interval_size = extended_interval[1] - extended_interval[0]
        assert interval_size + 100 == extended_interval_size

        # Extract the sequence in the interval. We'll modify this extracted
        # sequence only, and then put all the parts back together again
        # once we're done.
        before_interval_seq = genome_record.seq[:extended_interval[0]]
        orig_interval_seq = genome_record.seq[extended_interval[0]:
                extended_interval[1]]
        after_interval_seq = genome_record.seq[extended_interval[1]:]

        # Make a copy to manipulate so that we can make a comparison afterward.
        interval_seq = orig_interval_seq[:]
        assert extended_interval_size == len(interval_seq), (
                "Evaluating interval %s: \n %d != %d" % (
                        str(interval),
                        extended_interval_size,
                        len(interval_seq)))

        # Perform the fix. March through each position, updating the interval
        # seq as necessary until that position is above the threshold.
        for interval_pos in range(extended_interval_size):
            # We only care about increasing GC of windows centered at positions
            # inside the original interval, not in the extended interval.
            if (interval_pos < half_window or
                    interval_pos >= extended_interval_size - half_window):
                continue
            # Calculate window coordinates in frame of interval_seq.
            window_start = interval_pos - half_window
            window_end = interval_pos + half_window
            window_seq = interval_seq[window_start:window_end]

            # Figure out which positions we can modify.
            pos_to_modify_list = range(window_start, window_end)
            pos_to_modify_list = [p for p in pos_to_modify_list if not
                    _interval_pos_to_global_pos(p, extended_interval) in
                    pos_blacklist_set]

            # Repeat while we have positions to modify, or until we achieve
            # the desired GC content.
            gc_content = GC(window_seq) / 100
            while (len(pos_to_modify_list) and gc_content <
                    gc_content_constraint_obj.local_window_lower_bound):
                # Choose the next position randomly to avoid bias.
                pos_to_modify_idx = random.randint(
                        0, len(pos_to_modify_list) - 1)
                pos_to_modify = pos_to_modify_list.pop(pos_to_modify_idx)
                current = interval_seq[pos_to_modify]
                if current.upper() in GC_bases:
                    continue
                new = GC_TRANSITION_TABLE[current]
                interval_seq = (interval_seq[:pos_to_modify] + new +
                        interval_seq[pos_to_modify + 1:])
                window_seq = interval_seq[window_start:window_end]
                gc_content = GC(window_seq) / 100

        # Put it back together.
        assert len(orig_interval_seq) == len(interval_seq), (
            "len before: %d | len after: %d" % (
                    len(orig_interval_seq), len(interval_seq)))
        genome_record.seq = (before_interval_seq + interval_seq +
                after_interval_seq)

        # Post-completion checks.
        assert len_genome_before_fix == len(genome_record), (
                "len before: %d | len after: %d" % (
                        len_genome_before_fix, len(genome_record)))


def _interval_pos_to_global_pos(interval_pos, interval):
    return interval[0] + interval_pos


def _calculate_feature_annotation_shadow(genome_record,
        changes_ok_in_feature_types=[]):
    """Returns set of positions that are covered by features.

    This is useful for cases where we want to modify positions that are not
    in the shadow, like in intergenic GC content fixing.
    """
    pos_blacklist_set = set()

    def _update_pos_blacklist(pos_range):
        """Helper function that updates the blacklist.
        """
        for pos in pos_range:
            pos_blacklist_set.add(pos)

    filtered_features = [f for f in genome_record.features if f.type not
            in changes_ok_in_feature_types]

    for feature in filtered_features:
        start = feature.location.start
        end = feature.location.end
        _update_pos_blacklist(range(start, end))

        if feature.type == 'CDS':
            if feature.strand == 1:
                pos_range = range(start - 20, end)
            else:
                pos_range = range(end, end + 20)
            _update_pos_blacklist(pos_range)

    return pos_blacklist_set
