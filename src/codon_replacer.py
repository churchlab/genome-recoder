"""
Objects for replacing synonymous codons.

Our rules/heuristics for codon replacement depend on distributions
as well as adherance to FeatureProfilers.
"""

from Bio.Alphabet import generic_dna
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature

from annotation_feature_types import FORBIDDEN_CODON_FEATURE_TYPE
import biopython_util
from biopython_util import get_codon_feature_location
from biopython_util import translate_custom
from refactor_config import SELENOCYSTEINE_CODON


class BaseCodonReplacer:
    """Base class for objects that handle substituting synonymous codons across
    a genome feature (e.g. CDS), while maintaining heuristic profile
    constraints.

    Remember, our overall goal is to remove codons entirely from a genome, and
    to do so with the least impact on heuristic profiles that we wish to
    maintain (e.g. GC content, secondary structure, codon rarity, etc.).
    """
    codons_to_remove = None
    codon_usage_memex = None
    profiler_list = None


    def __init__(
            self,
            codons_to_remove,
            codon_usage_memex,
            profiler_list=[]):
        """Default constructor.

        Args:
            codons_to_remove: The list of codons we are removing from the
                Genome.
            codon_usage_memex: Object for looking up codon_usage.
            profiler_list: List of Profiler objects that we will use for
                enforcing constraints during synonymous codon substitution.
        """
        self.codons_to_remove = codons_to_remove
        self.codon_usage_memex = codon_usage_memex
        self.profiler_list = profiler_list


    def set_profiler_list(self, profiler_list):
        """Sets the profiler list.

        A CodonReplacer object can generally be re-used after just swapping out
        the profiler list.
        """
        self.profiler_list = profiler_list


    def replace_codons_in_feature(self, feature_id, seq_record):
        """Replaces the codons in the given feature relative to the given
        SeqRecord. Returns an object containing the new sequence. Does not
        mutate seq_record.

        Subclasses should override this.
        """
        feature = biopython_util.get_feature_by_id(seq_record, feature_id)
        feature_seq = seq_record.seq[
                feature.location.start:feature.location.end]

        return {
                'feature_id': feature_id,
                'orig_feature_seq': feature_seq,
                'new_feature_seq': feature_seq,
        }


class SimpleCodonReplacer(BaseCodonReplacer):
    """Replaces codons with any synonymous codon that we are not removing.

    We don't actually use this implementation, but it serves as a simpler
    demonstration on the idea of removing codons which should be understood
    before looking at the GraphSearchCodonReplacer.
    """

    def replace_codons_in_feature(self, feature_id, seq_record):
        feature = biopython_util.get_feature_by_id(seq_record, feature_id)
        feature_seq = seq_record.seq[
                feature.location.start:feature.location.end]

        new_feature_seq = ''
        for codon_start in range(0, len(feature_seq), 3):
            old_codon = str(feature_seq[codon_start:codon_start + 3])
            if old_codon not in self.codons_to_remove:
                new_feature_seq += old_codon
            else:
                new_codon = self.get_codon_replacement(old_codon)
                new_feature_seq += new_codon

        return {
                'is_success': True,
                'feature_id': feature_id,
                'orig_feature_seq': feature_seq,
                'new_feature_seq': Seq(new_feature_seq, generic_dna),
        }


    def get_codon_replacement(self, old_codon):
        """Returns a single codon to replace the given one.
        """
        synonymous_codons = self.codon_usage_memex.get_synonymous_codons(
                old_codon)
        for codon_name in synonymous_codons:
            if not codon_name in self.codons_to_remove:
                return codon_name
        raise RuntimeError("No valid codons to swap with %s" % old_codon)



class GraphSearchCodonReplacer(BaseCodonReplacer):
    """Naive implementation that does a full graph search to optimize the
    replacement scores for the profilers.

    We might want to add some basic pruning.
    """

    def replace_codons_in_feature(
            self,
            feature_id,
            seq_record,
            start_codon_index=0,
            last_codon_index=None,
            avoid_codons_in_positions={}):
        """Override default implementation to perform search over the feature
        profiles.

        NOTE: Let's avoid making full copies of the SeqRecord and insead just
        copy the feature sequence for now, then create the new
        SeqRecord at the end.

        Args:
            feature_id: Id of the feature to fix.
            seq_record: Source SeqRecord object. Not mutated.
            start_codon_index: Position to start mutating at. Useful for
                partial replacaement such as when resolving homology issues.
            avoid_codons_in_positions: Map from position to codon to use as
                last resort. Useful for resolving homology issues.

        Returns:
            An object with keys:
                * feature_id
                * is_success
                * exception_string
                * orig_feature_seq
                * new_feature_seq
                * profiler_stats
                * new_feature_list: List of features created during
                    refactoring. For example, we create a feature for
                    every forbidden codon created.
        """
        # Extract the start feature.
        feature = biopython_util.get_feature_by_id(
                seq_record, feature_id)

        # Extract the original feature sequence, using the correct polarity
        # depending on the strand.
        orig_feature_seq = seq_record.seq[
                feature.location.start:feature.location.end]
        if feature.strand == -1:
            orig_feature_seq = orig_feature_seq.reverse_complement()

        total_codons = len(orig_feature_seq) / 3
        if not last_codon_index:
            last_codon_index = total_codons - 1

        # Find the new feature sequence, watching for exceptions.
        is_success = True
        exception_string = None
        try:
            dfs_result = self._dfs(
                    feature,
                    orig_feature_seq,
                    start_codon_index,
                    last_codon_index,
                    avoid_codons_in_positions)
            new_feature_seq = dfs_result['feature_seq']
            new_feature_list = dfs_result['new_feature_list']
        except AssertionError as e:
            is_success = False
            exception_string = str(e)
            new_feature_seq = None
            new_feature_list = []

        # Save some metadata stats from this calculation.
        profiler_stats = {}
        for profiler in self.profiler_list:
            profiler_stats[profiler.__class__.get_name()] = {
                'error_tolerance': profiler.error_tolerance
            }

        return {
                'feature_id': feature_id,
                'is_success': is_success,
                'exception_string': exception_string,
                'orig_feature_seq': orig_feature_seq,
                'new_feature_seq': new_feature_seq,
                'profiler_stats': profiler_stats,
                'new_feature_list': new_feature_list
        }


    def _dfs(
            self,
            orig_feature,
            orig_feature_seq,
            start_codon_index,
            last_codon_index,
            avoid_codons_in_positions):
        """Perform depth-first search to get to the solution that works.

        The strategy to is to maintain a queue ordered by nodes to explore.
        When we extend a node, we add the resulting nodes to the front of the
        queue. If the profile score ever becomes sufficiently bad for a
        particular node we throw it on the end of the queue.
        """
        # NOTE: It seems that a 'while True' is okay here, since every time we
        # we get stuck we increment the profile error tolerance and repeat>
        # (Eventually the profile error_tolerances will be silly high so that
        # any solution should work, but we should probably think about placing
        # a limit on this somehow.)
        while True:
            queue = []

            # Each node is an object with a sequence and a pointer to
            # the next codon that needs to be changed.
            start_node = {
                    'feature_seq': orig_feature_seq,
                    'next_codon_index': start_codon_index,
                    'profiler_score_obj_list': [0 for p in self.profiler_list],
                    'new_feature_list': []
            }
            queue.append(start_node)

            total_codons = len(orig_feature_seq) / 3

            # Because some sequences may be particularly tricky with regard to
            # adhering to the error_tolerance of a feature_profile, we want to
            # be able to give these troublesome sequences some room to breath
            # by incrementing the error_tolerances for those profiles to some
            # extent. Thus, here we need to somehow detect that the search is
            # getting stuck. A first attempt strategy for accomplishing this is as
            # follows:
            # We maintain two state variables:
            #     * max_index_observed: Keep track of the furthest codon we've
            #           tried to fix so far.
            #     * iterations_at_max: Track how long it's been since we've seen
            #           a new max. If this goes beyond some threshold, interrupt
            #           the loop, increment the profile error_tolerances,
            #           and start searching from the beginning again.

            ITERATION_LIMIT = 200

            max_index_observed = 0
            iterations_at_max = 0

            # Pop the first item off the queue and extend it.
            while queue:
                current_node = queue.pop(0)

                # Logic for tracking whether we are probably stuck.
                current_index = current_node['next_codon_index']
                if current_index > max_index_observed:
                    max_index_observed = current_index
                    iterations_at_max = 0
                else:
                    iterations_at_max += 1

                if iterations_at_max >= ITERATION_LIMIT:
                    # print '...Detected that search got stuck'
                    break

                # print 'Replacing codon %d of %d' % (
                #         current_index + 1, total_codons)

                # If this node is done, we've reached the end, return it.
                if current_index > last_codon_index:
                    self._final_node_validation(
                            current_node, orig_feature, orig_feature_seq)
                    return current_node

                self._validate_node_codon(current_node, total_codons)

                # Extend the current node.
                extended_nodes = self._extend_node(orig_feature, current_node,
                        total_codons, avoid_codons_in_positions)

                # Since this is a depth-first search, add the results to the
                # front of the queue.
                queue = extended_nodes + queue

            # If we're here, that means we couldn't find a solution with the
            # present feature_profile error tolerances, so we:
            #     * Identify the worst one and increment its error tolerance
            #     * Reset all profile error counts
            # And then return to the beginning of the while loop.
            most_failed_profile = self.profiler_list[0]
            most_failed_count = most_failed_profile.failed_error_treshold_count
            for profile in self.profiler_list[1:]:
                current_failed_count = profile.failed_error_treshold_count
                if current_failed_count > most_failed_count:
                    most_failed_profile = profile
                    most_failed_count = current_failed_count
            most_failed_profile.increment_error_tolerance()

            # Reset the failed error counts.
            for profile in self.profiler_list:
                profile.reset_failed_error_threshold_count()


    def _validate_node_codon(self, node, total_codons):
        """Makes sure that the codon of the current node is valid.

        Basically, we're just checking that we don't see a stop codon
        in a position that is not the last one, and if it is the last
        position, then the codon must be a stop codon.

        Raises:
            AssertionError if this is an invalid codon for this position.
        """
        codon_start_pos = node['next_codon_index'] * 3
        current_codon = str(
                node['feature_seq'][codon_start_pos:codon_start_pos + 3])

        # If this is the last one, it must be a stop codon, otherwise
        # it must not be one.
        if node['next_codon_index'] == total_codons - 1:
            assert self.codon_usage_memex.is_stop_codon(current_codon),\
                    'Last codon is not stop codon.'
        elif SELENOCYSTEINE_CODON == current_codon:
            # TGA occurring mid-feature, codes for Selenocysteine.
            return True
        else:
            assert not self.codon_usage_memex.is_stop_codon(current_codon),\
                    'Stop codon observed at codon index %d' % (
                            node['next_codon_index'],)
        return True


    def _extend_node(self, orig_feature, node, total_codons,
            avoid_codons_in_positions):
        """Pick synonymous codons that are not forbidden for the next
        codon.  Perform profile-based scoring.  Order them according to
        best score so far.

        Args:
            orig_feature: The original feature we are making changes
                inside of.
            node: Object with keys:
                * feature_seq
                * next_codon_index
                * profiler_score_obj_list
            total_codons: Number of codons in this feature.
            avoid_codons_in_positions: Map from position to codon to use as
                last resort. Useful for resolving homology issues.

        Returns:
            List of node objects.
        """
        codon_index = node['next_codon_index']
        codon_start_pos = codon_index * 3
        current_codon = str(
                node['feature_seq'][codon_start_pos:codon_start_pos + 3])

        # Don't replace, just move on, if:
        #    * We are looking at the start codon (we don't replace start codons).
        #    * Codon is TGA and we are not at the end (Selenocysteine).
        #    * Non-forbidden codon, unless we're trying to explicitly
        #      avoid this codon at this position by indicating it in the
        #      avoid_codons_in_positions map.
        if ((codon_index == 0) or
                (current_codon == SELENOCYSTEINE_CODON and
                        codon_index != total_codons - 1) or
                (current_codon not in self.codons_to_remove and
                        codon_index not in avoid_codons_in_positions)):
            new_node = {
                    'feature_seq': node['feature_seq'][:],
                    'next_codon_index': node['next_codon_index'] + 1,
                    'profiler_score_obj_list': node['profiler_score_obj_list'],
                    'new_feature_list': node['new_feature_list']
            }
            return [new_node]

        if codon_start_pos == 0:
            synonymous_codons = self.codon_usage_memex.get_start_codon_list()
        else:
            synonymous_codons = self.codon_usage_memex.get_synonymous_codons(
                current_codon)
        assert current_codon in synonymous_codons

        new_nodes = []
        for codon_name in synonymous_codons:
            if codon_name in self.codons_to_remove:
                continue

            # Avoid codon if there are other options.
            if (codon_index in avoid_codons_in_positions):
                try_avoid_codon = avoid_codons_in_positions[codon_index]
                if (set(synonymous_codons) -
                        set(self.codons_to_remove) -
                        set([try_avoid_codon])):
                    if (codon_name == try_avoid_codon):
                        continue

            # At this point, it's possible that we are just swapping in the
            # same codon if we could not avoid it (e.g. it was the only
            # option).

            new_feature_seq = Seq(
                    str(node['feature_seq'][:codon_start_pos]) +
                            codon_name +
                            str(node['feature_seq'][codon_start_pos + 3:]),
                    generic_dna
            )

            # Calculate errors relative to original profiles. All profiles
            # must pass in order for the node to be eligible for further
            # expansion.
            profiler_score_obj_list = []
            all_passed = True
            for profiler in self.profiler_list:
                score_obj = profiler.calc_feature_score_relative_to_profile(
                        node['next_codon_index'], new_feature_seq)
                if not score_obj['passes_error_threshold']:
                    all_passed = False
                    break
                profiler_score_obj_list.append(score_obj)

            if not all_passed:
                continue

            # If we're changing the codon, add a feature annotating the change.
            new_feature_list = node['new_feature_list'][:]
            if current_codon != codon_name:
                new_feature_id = (str(orig_feature.id) + '_' +
                        FORBIDDEN_CODON_FEATURE_TYPE +
                        '_' + str(codon_index) +
                        '_' + current_codon +
                        '_' + codon_name)
                new_feature_location_start = (orig_feature.location.start +
                        codon_start_pos)
                new_feature_location = get_codon_feature_location(orig_feature,
                        codon_index)
                new_feature = SeqFeature(
                        location=new_feature_location,
                        type=FORBIDDEN_CODON_FEATURE_TYPE,
                        strand=orig_feature.strand,
                        id=new_feature_id)

                new_feature.qualifiers['ref'] = current_codon

                new_feature.qualifiers['alt'] = codon_name

                new_feature.qualifiers['delta_GC'] = str(
                        calc_delta_GC_count(current_codon, codon_name))

                # Add delta secondary structure free energy.
                for score_obj in profiler_score_obj_list:
                    if score_obj['profile_name'] == 'ss_profile':
                        new_feature.qualifiers['delta_ss'] = str(
                                score_obj['error_score'])
                        break

                new_feature_list += [new_feature]

            # Package up the data for the new node in the search graph.
            new_node = {
                    'feature_seq': new_feature_seq,
                    'next_codon_index': node['next_codon_index'] + 1,
                    'profiler_score_obj_list': profiler_score_obj_list,
                    'profiler_ranks': [], # Used for sorting below.
                    'new_feature_list': new_feature_list
            }
            new_nodes.append(new_node)

        # At this point we might want to sort the results according to
        # some reasonable measure so that we increase our chances of going down
        # good paths. Remember, since we are performing a depth-first search,
        # this node.

        # Current strategy: Come up with a ranking of each node relative
        # to the others in terms of scores relative to the profiles. Then sort
        # according to a weighted function of "importance" of that feature.
        for profiler_index in range(len(self.profiler_list)):
            profiler = self.profiler_list[profiler_index]
            nodes_ordered_for_this_profile = sorted(
                    new_nodes,
                    key=lambda x: x['profiler_score_obj_list'][profiler_index]
            )
            for node_rank in range(len(nodes_ordered_for_this_profile)):
                node = nodes_ordered_for_this_profile[node_rank]
                node['profiler_ranks'].append(node_rank + 1)
        for node in new_nodes:
            node['weighted_rank'] = sum(node['profiler_ranks'])


        # Sort the nodes, lower weighted_rank first.
        sorted_new_nodes = sorted(
                new_nodes,
                key=lambda x: x['weighted_rank']
        )

        return sorted_new_nodes


    def _final_node_validation(self, node, orig_feature, orig_feature_seq):
        """Validate the resulting sequence.

        Make sure:
            * Same translation (amino acid sequence preserved).

        Raises:
            AssertionError if any of the validations fail.
        """
        try:
            original_translation_str = str(translate_custom(orig_feature_seq))
            new_translation_str = str(translate_custom(node['feature_seq']))
        except TranslationError as e:
            copyargs = list(e.args)
            copyargs[0] += '- Feature id: %s' % str(orig_feature.id)
            e.args = tuple(copyargs)
            raise e
        assert original_translation_str == new_translation_str,\
                "Original translation %s does not match new translation %s" % (
                        original_translation_str, new_translation_str)

GC_SET = set('GC')

def calc_delta_GC_count(original_seq, new_seq):
    """Returns the difference in GC count between two sequences.

    For simplicity, assume the two sequences are of the same length.
    """
    def _count_gc(seq):
        """Returns the number of G's or C's in the sequence.
        """
        return reduce(
                lambda count, base: count + 1 if base in GC_SET else count,
                seq, 0)

    return _count_gc(new_seq) - _count_gc(original_seq)
