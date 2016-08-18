"""
Methods for separating out overlapping genes.
"""

import pdb

import csv
import copy
import os
import pickle
import re

from Bio.Seq import reverse_complement
from Bio.Seq import translate
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.Data.CodonTable import TranslationError

from annotation_feature_types import FIX_OVERLAP_SIZE_4
from annotation_feature_types import FIX_OVERLAP_SYNONYMOUS
from biopython_util import add_feature_to_seq_record
from biopython_util import COMMONLY_IGNORED_FEATURE_TYPES
from biopython_util import get_feature_by_id
from biopython_util import get_feature_gene
from biopython_util import get_region_codon_indeces_in_feature
from biopython_util import insert_sequence
from biopython_util import update_feature_seq
from biopython_util import InsertType
from biopython_util import shift_seq_feature_right_with_head_copy
from biopython_util import replace_feature
from biopython_util import swap_feature_codon_at_position
from biopython_util import translate_custom
from conflicting_pair_finder import find_all_overlaps
from conflicting_pair_common import does_feature_have_codon_in_list_region
from conflicting_pair_common import does_feature_have_forbidden_codon_in_region
from conflicting_pair_common import does_upstream_feature_overlap_downstream_feature
from conflicting_pair_common import identify_conflict_region_of_interest
from conflicting_pair_common import is_feature_in_region
from conflicting_pair_common import RBS_BUFFER_SIZE
from refactor_config import CACHE_DIR
from tool_wrappers.rbs_calc_util import get_mds42_rbs_strength_profile
from tool_wrappers.rbs_calc_util import is_rbs_strength_conserved

FIX_OVERLAPS_CACHE_LOCATION = os.path.join(
        CACHE_DIR, 'mds42_fix_overlaps.cache')

AGN_CODONS = set(['AGA', 'AGC', 'AGG', 'AGT'])

class ConflictingPairFixer(object):
    """Object that encapsulates logic.

    NOTE: Currently in transition to this cleaner class-based design.
    """
    def __init__(self,
            refactor_context,
            cache=False,
            include_close_features=False,
            single_iteration=False,
            force_separate_AGN=False,
            agn_separation_data_file=None,
            debug_report_output='AGN_debug.csv'):
        # Basic validation.
        if force_separate_AGN:
            assert agn_separation_data_file, "If separating AGN, need data."

        # Assign attributes
        self.genome_record = refactor_context.get_genome_record()
        self.forbidden_codons = refactor_context.get_forbidden_codon_set()
        self.original_codon_usage_memex = (
                refactor_context.get_original_codon_usage_memex())
        self.cache = cache
        self.include_close_features = include_close_features
        self.single_iteration = single_iteration
        self.force_separate_AGN = force_separate_AGN
        self.debug_report_output = debug_report_output

        # Calculate others.
        self.original_rbs_strength_profile = get_mds42_rbs_strength_profile()

        # Maybe prepare AGN separation data.
        if agn_separation_data_file:
            self.agn_separation_data = {}
            with open(agn_separation_data_file) as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    assert not row['ID'] in self.agn_separation_data
                    self.agn_separation_data[row['ID']] = row

    def fix_overlaps(self):
        """Fix any overlapping features in essential_features, propagating
        sequence and position changes throughout the genome_record SeqRecord
        object and all of its features (not just the ones in essential_features).

        Algorithm overview:
            Initial check for overlapping pairs.
            While there is at least one overlapping pair remaining:
                Fix the particular overlap.
                Check for pairs again.

        Args:
            genome_record: The SeqRecord object that contains all the data.
                The result SeqRecord returned will be a derivative of this.
            cache: If True, look for existing cached results, only doing
                calculation if not found, and caching the results for the
                next run.

        Returns:
            An updated SeqRecord object with all overlaps resolved.
        """
        print 'Fixing overlaps ...'

        # NOTE: Temporary manual assignment of names during refactor.
        forbidden_codons = self.forbidden_codons
        original_codon_usage_memex = self.original_codon_usage_memex
        cache = self.cache
        include_close_features = self.include_close_features
        single_iteration = self.single_iteration

        # First try to recover cached results:
        if cache:
            try:
                with open(FIX_OVERLAPS_CACHE_LOCATION) as fh:
                    fixed_overlaps_genome_record = pickle.load(fh)
                    print ('...Using cached fixed overlaps from ' +
                            FIX_OVERLAPS_CACHE_LOCATION)
                    return fixed_overlaps_genome_record
            except:
                # Cache doesn't exist. Add more precise error catching if desired.
                pass


        ### Some initial validation and sanity checks.

        # Make sure ids are unique.
        error_msg = (
                "overlaps.fix_overlaps() requires all SeqFeatures to have a "
                "unique id")
        id_set = set(map(lambda f: f.id, self.genome_record.features))
        assert len(id_set) == len(self.genome_record.features), error_msg

        def _find_overlaps(current_genome_record, cache=False):
            """Internal method to keep the code in one place.
            Though we only cache on the first call.
            """
            overlaps = find_all_overlaps(
                    current_genome_record,
                    forbidden_codons,
                    ignore_feature_type=COMMONLY_IGNORED_FEATURE_TYPES,
                    full_analysis=False,
                    cache=cache,
                    include_close_features=include_close_features)
            return overlaps

        overlaps = _find_overlaps(self.genome_record, cache=cache)

        # Set of ids for overlap pairs.
        fixed_overlaps = set()

        debug_report_rows = []

        # Keep track of iterations for debugging.
        iteration = 1
        while overlaps:
            num_overlaps = len(overlaps)
            print 'Fixing overlaps iteration ' + str(iteration)
            for overlap_index, conflict_param_obj in enumerate(overlaps):
                print 'Fixing overlap %d of  %d' % (overlap_index + 1, num_overlaps)
                upstream_feature_id = conflict_param_obj['upstream_feature'].id
                downstream_feature_id = conflict_param_obj['downstream_feature'].id

                upstream_feature = conflict_param_obj['upstream_feature']
                u_seq = upstream_feature.extract(self.genome_record.seq)
                print 'upstream', upstream_feature_id
                print u_seq[:20], '...', u_seq[-20:]

                downstream_feature = conflict_param_obj['downstream_feature']
                print 'downstream', downstream_feature_id
                d_seq = downstream_feature.extract(self.genome_record.seq)
                print d_seq[:20], '...', d_seq[-20:]

                if include_close_features:
                    fix_conservative_result = (
                            self.fix_conflicting_pair_conservative(
                                    conflict_param_obj,
                                    self.genome_record,
                                    forbidden_codons,
                                    original_codon_usage_memex
                            )
                    )
                    is_fixed = fix_conservative_result['is_success']
                    if 'debug_report_row' in fix_conservative_result:
                        debug_report_rows.append(
                                fix_conservative_result['debug_report_row'])
                    if is_fixed:
                        # NOTE: This is a NOOP but safer to leave it
                        # in case we do away with mutating methods.
                        self.genome_record = fix_conservative_result[
                                'updated_seq_record']
                else:
                    is_fixed = fix_overlap_pair(
                            upstream_feature_id,
                            downstream_feature_id,
                            self.genome_record)
                assert is_fixed, "Could not fix pair %s and %s" % (
                        upstream_feature_id, downstream_feature_id)

                # Track overlaps that were fixed.
                fixed_overlaps.add(self._make_overlap_id(conflict_param_obj))

            # Check again in case new overlaps were introduced somehow.
            if single_iteration:
                overlaps = []
                DEBUG_REPORT_FIELD_NAMES = [
                        'feature_id',
                        'feature_gene',
                        'strand',
                        'LP',
                        'RP',
                        'seq',
                ]
                if len(debug_report_rows) > 0:
                    with open(self.debug_report_output, 'w') as fh:
                        writer = csv.DictWriter(fh, DEBUG_REPORT_FIELD_NAMES)
                        for row in debug_report_rows:
                            writer.writerow(row)
            else:
                # NOTE: Follow up check should not be cached as fixing it may
                # affect it.
                overlaps = _find_overlaps(self.genome_record, cache=False)
                overlaps = filter(
                        lambda pair_obj:
                                (not self._make_overlap_id(pair_obj) in
                                        fixed_overlaps),
                        overlaps
                )

            iteration += 1

        # Maybe cache the results.
        if cache:
            if not os.path.exists(CACHE_DIR):
                os.mkdir(CACHE_DIR)
            with open(FIX_OVERLAPS_CACHE_LOCATION, 'w') as fh:
                pickle.dump(self.genome_record, fh)

        # The resulting SeqRecord should no longer have overlaps among
        # the features with ids in essential_feature_ids.
        return self.genome_record


    def _make_overlap_id(self, pair_obj):
        """Helper methods that makes a pair id from an objecting representing
        a pair.
        """
        upstream_feature_id = str(pair_obj['upstream_feature'].id)
        downstream_feature_id = str(pair_obj['downstream_feature'].id)
        return upstream_feature_id + '__' + downstream_feature_id


    def fix_conflicting_pair_conservative(
            self,
            conflict_param_obj,
            seq_record,
            forbidden_codons,
            codon_usage_memex):
        """Finds the most conservative way to resolve a conflict between two
        features.

        NOTE: This method mutates seq_record, specifically the underlying
        sequence and any feature positions as necessary.

        Algorithm overview:
            * If either of the affected regions is not a coding feature:
                - It's not possible to resolve it using synonymous swaps,
                  continue to forceful separation step.
            * If there are no forbidden codons in the affected regions:
                - Nothing to do, mark the pair as resolved.
            * If there are forbidden codons:
                - Do an exhaustive search over the affected region
                  and try to find a path of synonymous codon
                  substitutions that don't require physical separation.
                  If success, perform the change, and mark any changed
                  codons as "fixed" so that they are not changed in the
                  second half of the overall algorithm where we do the bulk
                  forbidden codon removal.
                - Otherwise, we need to separate:
                    * If overlap < 4 bp, find minimum amount to copy that
                      resolves any issues, and lock affected RBS regions
                      so that they are not changed.
                    * Otherwise, need to copy overlap + 15 bp upstream
                      of ATG, and to help prevent snap-back:
                        - muddle old start codon in upstream gene
                        - muddle bases in copied region that are not part
                          of RBS

        Returns:
            Object with keys:
                * is_success: Boolean indicating whether the conflict is
                        resolved.
                * updated_seq_record: If is_success is True.

        Raises:
            AssertionError: If there is no conflict to fix.
        """
        fix_result = {
                'is_success': False
        }

        # Controls extent to which we try conservative fix.
        MAX_OVERLAP_SIZE_TO_FIX = 100

        # Parse params.
        upstream_feature_id = conflict_param_obj['upstream_feature'].id
        downstream_feature_id = conflict_param_obj['downstream_feature'].id
        conflict_type = conflict_param_obj['conflict_type']

        upstream_feature = get_feature_by_id(
                seq_record, upstream_feature_id)
        assert upstream_feature
        print '...UPSTREAM_FEATURE', upstream_feature.id, upstream_feature.location

        downstream_feature = get_feature_by_id(
                seq_record, downstream_feature_id)
        assert downstream_feature
        print '...DOWNSTREAM_FEATURE', downstream_feature.id, downstream_feature.location

        # If both features are not coding features, we can't do synonymous
        # swaps (might break RNA) and so proceed directly to hard overlap fix.
        if upstream_feature.type != 'CDS' or downstream_feature.type != 'CDS':
            is_success = fix_overlap_pair(
                    upstream_feature_id, downstream_feature_id, seq_record)
            fix_result['is_success'] = is_success
            if fix_result['is_success']:
                fix_result['updated_seq_record'] = seq_record
            return fix_result

        # First, identify the relevant region of the genome to look at. That is,
        # we'll look in the affected features in this region.
        region_of_interest = identify_conflict_region_of_interest(
                upstream_feature, downstream_feature)
        if not region_of_interest:
            # Nothing to do. Mark as fixed.
            fix_result['is_success'] = True
            fix_result['updated_seq_record'] = seq_record
            return fix_result

        # # Repeat the forbidden codon check performed while finding overlaps.
        # # This is a cheap check, so leaving it in to catch any bugs that might
        # # creep in.
        # upstream_has_forbidden = (
        #         upstream_feature.type == 'CDS' and
        #         is_feature_in_region(upstream_feature, region_of_interest) and
        #         does_feature_have_forbidden_codon_in_region(
        #                 upstream_feature, region_of_interest, seq_record,
        #                 forbidden_codons)
        # )
        # downstream_has_forbidden = (
        #         downstream_feature.type == 'CDS' and
        #         is_feature_in_region(downstream_feature, region_of_interest) and
        #         does_feature_have_forbidden_codon_in_region(
        #                 downstream_feature, region_of_interest, seq_record,
        #                 forbidden_codons)
        # )
        # if not (upstream_has_forbidden or downstream_has_forbidden):
        #     assert False, ("Finding overlaps should not report an " +
        #             "overlapping region with no forbidden codons.")

        if self.force_separate_AGN:
            # If at least one AGN codon in region of interest, then perform
            # hard fix.
            has_relevant_AGN = False
            debug_report_row = None
            if upstream_feature.strand == 1 and downstream_feature.strand == 1:
                codon_index = does_feature_have_codon_in_list_region(
                        upstream_feature, region_of_interest, seq_record,
                        AGN_CODONS, return_codon_index=True)
                if codon_index:
                    if upstream_feature.id in self.agn_separation_data:
                        separation_data = self.agn_separation_data[
                                upstream_feature.id]
                        if separation_data['Separate'] == '1':
                            has_relevant_AGN = True
                        # Manually fix the codon.
                        if 'Codon' in separation_data:
                            new_codon = separation_data['Codon'].upper()
                            if len(new_codon) == 3:
                                upstream_feature_seq = str(upstream_feature.extract(
                                        seq_record.seq))
                                original_codon = upstream_feature_seq[
                                        codon_index * 3:codon_index * 3 + 3]
                                assert original_codon in AGN_CODONS
                                seq_record = swap_feature_codon_at_position(
                                        seq_record, upstream_feature.id,
                                        codon_index * 3, original_codon, new_codon)

                #     codon_index -= 2
                #     seq_string = ''
                #     genome_pos = (upstream_feature.location.start +
                #             codon_index * 3)
                #     LP = genome_pos
                #     while genome_pos < downstream_feature.location.start:
                #         if genome_pos + 3 >= downstream_feature.location.start:
                #             # We're going to bleed over the start of the
                #             # next feature so write it all together.
                #             remaining = (seq_record.seq[genome_pos:
                #                     downstream_feature.location.start]).lower()
                #             next_start_codon = str(seq_record.seq[
                #                     downstream_feature.location.start:
                #                     downstream_feature.location.start + 3])
                #             seq_string += remaining + next_start_codon
                #             genome_pos = downstream_feature.location.start + 3
                #             break
                #         codon = str(seq_record.seq[genome_pos:
                #                 min(upstream_feature.location.end, genome_pos + 3)
                #         ])
                #         if not codon in self.forbidden_codons:
                #             codon = codon.lower()
                #         gap_length = 3 - len(codon)
                #         if gap_length == 0:
                #             seq_string += codon + ' '
                #         else:
                #             seq_string += codon + ''.join(['*'] * gap_length)
                #         genome_pos += 3
                #     RP = genome_pos
                #     if len(seq_string) > 0:
                #         return {
                #             'is_success': True,
                #             'debug_report_row': {
                #                 'feature_id': upstream_feature.id,
                #                 'feature_gene': get_feature_gene(
                #                         upstream_feature),
                #                 'strand': 1,
                #                 'LP': LP,
                #                 'RP': RP,
                #                 'seq': seq_string
                #             }
                #         } # DO NOT SUBMIT
                #     else:
                #         return {'is_success': True}

            elif upstream_feature.strand == -1 and downstream_feature.strand == -1:
                codon_index = does_feature_have_codon_in_list_region(
                        downstream_feature, region_of_interest, seq_record,
                        AGN_CODONS, return_codon_index=True)
                if codon_index:
                    if downstream_feature.id in self.agn_separation_data:
                        separation_data = self.agn_separation_data[
                                downstream_feature.id]
                        if separation_data['Separate'] == '1':
                            has_relevant_AGN = True
                        # Manually fix the codon.
                        if 'Codon' in separation_data:
                            new_codon = separation_data['Codon'].upper()
                            if len(new_codon) == 3:
                                downstream_feature_seq = str(
                                        downstream_feature.extract(seq_record.seq))
                                original_codon = downstream_feature_seq[
                                        codon_index * 3:codon_index * 3 + 3]
                                assert original_codon in AGN_CODONS, (
                                        original_codon + str(separation_data))
                                new_codon = separation_data['Codon'].upper()
                                seq_record = swap_feature_codon_at_position(
                                        seq_record, downstream_feature.id,
                                        codon_index * 3, original_codon, new_codon)

                    # codon_index -= 2
                    # downstream_feature_seq = downstream_feature.extract(
                    #         seq_record.seq)
                    # seq_string = ''
                    # genome_pos = (downstream_feature.location.end - 1 -
                    #         codon_index * 3)
                    # RP = genome_pos
                    # while genome_pos >= upstream_feature.location.end:
                    #     if genome_pos - 3 <= upstream_feature.location.end:
                    #         # We're going to bleed over the start of the
                    #         # next feature so write it all together.
                    #         remaining = str(reverse_complement(
                    #             seq_record.seq[
                    #                     upstream_feature.location.end:
                    #                     genome_pos + 1])).lower()
                    #         next_start_codon = str(reverse_complement(
                    #             seq_record.seq[
                    #                     upstream_feature.location.end - 3:
                    #                     upstream_feature.location.end]))
                    #         seq_string += remaining + next_start_codon
                    #         genome_pos = upstream_feature.location.end - 1 - 3
                    #         break
                    #     codon = str(downstream_feature_seq[codon_index * 3:
                    #             codon_index * 3 + 3])
                    #     if not codon in self.forbidden_codons:
                    #         codon = codon.lower()
                    #     gap_length = 3 - len(codon)
                    #     if gap_length == 0:
                    #         seq_string += codon + ' '
                    #     else:
                    #         seq_string += codon + ''.join(['*'] * gap_length)
                    #     codon_index += 1
                    #     genome_pos = (downstream_feature.location.end - 1 -
                    #         codon_index * 3)
                    # LP = genome_pos
                    # if len(seq_string) > 0:
                    #     return {
                    #         'is_success': True,
                    #         'debug_report_row': {
                    #             'feature_id': downstream_feature.id,
                    #             'feature_gene': get_feature_gene(
                    #                     downstream_feature),
                    #             'strand': -1,
                    #             'LP': LP,
                    #             'RP': RP,
                    #             'seq': seq_string
                    #         }
                    #     } # DO NOT SUBMIT
                    # else:
                    #     return {'is_success': True}

            elif upstream_feature.strand == -1 and downstream_feature.strand == 1:
                # Ignore this case for now.
                has_relevant_AGN = False
                # has_relevant_AGN = does_feature_have_codon_in_list_region(
                #         upstream_feature, region_of_interest, seq_record,
                #         AGN_CODONS)
                # has_relevant_AGN = (has_relevant_AGN or
                #         does_feature_have_codon_in_list_region(
                #                 downstream_feature, region_of_interest, seq_record,
                #                 AGN_CODONS))
            if has_relevant_AGN:
                print '...FORCING SPLIT DUE TO AGN IN RBS.'
                is_success = fix_overlap_pair(
                        upstream_feature_id, downstream_feature_id, seq_record)
                fix_result['is_success'] = is_success
                if fix_result['is_success']:
                    fix_result['updated_seq_record'] = seq_record
                return fix_result

        # Limit how big a region we try to tackle. If it's too big it becomes
        # computationally overwhelming.
        region_size = region_of_interest[1] - region_of_interest[0]
        if region_size > MAX_OVERLAP_SIZE_TO_FIX:
            is_success = fix_overlap_pair(
                    upstream_feature_id, downstream_feature_id, seq_record)
            fix_result['is_success'] = is_success
            if fix_result['is_success']:
                fix_result['updated_seq_record'] = seq_record
            return fix_result

        print '...REGION:', region_of_interest, region_size

        # Since we have forbidden codons, first we try to resolve them using
        # synonymous swaps only.
        synonymous_swap_result = (
            fix_conflict_using_synonymous_codon_resolution(
                    upstream_feature,
                    downstream_feature,
                    conflict_type,
                    region_of_interest,
                    seq_record,
                    forbidden_codons,
                    codon_usage_memex,
                    self.original_rbs_strength_profile
            ))
        if synonymous_swap_result['is_success']:
            # Add feature annotation for this resolution.
            feature_type = FIX_OVERLAP_SYNONYMOUS
            feature_id = ('%s_%s_%s' % (
                    feature_type,
                    get_feature_gene(upstream_feature),
                    get_feature_gene(downstream_feature)))
            feature_location = FeatureLocation(
                    region_of_interest[0], region_of_interest[1])
            feature = SeqFeature(
                    type=feature_type,
                    location=feature_location,
                    strand=1,
                    id=feature_id
            )
            add_feature_to_seq_record(
                    synonymous_swap_result['updated_seq_record'], feature)

            return synonymous_swap_result

        # Overlaps of size 4 that aren't resolved yet are a special case
        # which we treat as follows:
        if (conflict_param_obj['conflict_type'] == 'overlap' and
                conflict_param_obj['overlap_size'] == 4):
            fix_result = fix_same_dir_size_4_overlap(
                    upstream_feature,
                    downstream_feature,
                    conflict_param_obj['conflict_type'],
                    region_of_interest,
                    seq_record,
                    forbidden_codons,
                    codon_usage_memex,
                    self.original_rbs_strength_profile
            )
            if fix_result['is_success']:
                return fix_result

        # Otherwise, we need to pull them apart and copy the respective RBS.
        is_success = fix_overlap_pair(
                upstream_feature_id, downstream_feature_id, seq_record)
        fix_result['is_success'] = is_success
        if fix_result['is_success']:
            fix_result['updated_seq_record'] = seq_record
        return fix_result



###############################################################################
# Static methods.
###############################################################################

def fix_conflict_using_synonymous_codon_resolution(
        upstream_feature,
        downstream_feature,
        conflict_type,
        conflict_region,
        seq_record,
        forbidden_codons,
        codon_usage_memex,
        original_rbs_strength_profile):
    """Method that attempts to resolve a conflict between two features using
    synonymous codon swaps only.
    """
    print '...Attempt fix conflict using synonymous codon swaps...'

    # Generate all feature seq variants that remove forbidden codons.
    # First try swapping forbidden codons only to see if that works,
    # but then allow swapping any synonymous codons.
    for swap_forbidden_only in [True, False]:
        print '......swap forbidden only:', str(swap_forbidden_only)

        upstream_feature_seq_variants = _generate_permitted_feature_seq_variants(
            upstream_feature,
            downstream_feature,
            conflict_type,
            seq_record,
            conflict_region,
            forbidden_codons,
            codon_usage_memex,
            swap_forbidden_only=swap_forbidden_only)
        print '......upstream variants: %d' % len(upstream_feature_seq_variants)
        if not upstream_feature_seq_variants:
            continue

        # Make sure the set generated is unique.
        assert (len(upstream_feature_seq_variants) ==
                len(set(upstream_feature_seq_variants)))

        downstream_feature_seq_variants = _generate_permitted_feature_seq_variants(
            downstream_feature,
            upstream_feature,
            conflict_type,
            seq_record,
            conflict_region,
            forbidden_codons,
            codon_usage_memex,
            swap_forbidden_only=swap_forbidden_only)
        print '......downstream variants: %d' % len(downstream_feature_seq_variants)
        if not downstream_feature_seq_variants:
            continue

        # Make sure the set generated is unique.
        assert (len(downstream_feature_seq_variants) ==
                len(set(downstream_feature_seq_variants)))

        # Try to find a valid combination while extending the allowed order of
        # magnitude error.
        for rbs_allowed_error_order_of_magnitude in [1, 2, 3, 4]:
            print ('......Iterating over variants with allowed rbs error: %d' %
                    rbs_allowed_error_order_of_magnitude)

            # Find any valid combination.
            # TODO: Find the "best" valid combo, (and define "best").
            for upstream_variant in upstream_feature_seq_variants:
                for downstream_variant in downstream_feature_seq_variants:
                    is_valid_combo = is_valid_feature_seq_combo(
                        upstream_feature,
                        upstream_variant,
                        downstream_feature,
                        downstream_variant,
                        seq_record,
                        original_rbs_strength_profile,
                        rbs_allowed_error_order_of_magnitude
                    )
                    if is_valid_combo:
                        update_feature_seq(
                                seq_record,
                                upstream_feature.id,
                                upstream_variant)
                        update_feature_seq(
                                seq_record,
                                downstream_feature.id,
                                downstream_variant)

                        return {
                                'is_success': True,
                                'updated_seq_record': seq_record
                        }

    return {
            'is_success': False
    }


def is_valid_feature_seq_combo(
        upstream_feature,
        upstream_variant,
        downstream_feature,
        downstream_variant,
        seq_record,
        original_rbs_strength_profile,
        rbs_allowed_error_order_of_magnitude):
    """Tests whether the feature sequence reassignment is valid.

    This includes checking that the RBS expression strength is conserved,
    as calculated by the Salis RBS calculator.
    """
    # Get the forward strands to check overlay.
    if upstream_feature.strand == 1:
        upstream_feature_forward_seq = upstream_variant
    else:
        upstream_feature_forward_seq = reverse_complement(
                upstream_variant)

    if downstream_feature.strand == 1:
        downstream_feature_forward_seq = downstream_variant
    else:
        downstream_feature_forward_seq = reverse_complement(
                downstream_variant)

    # After checking overlay if necessary, we'll compuate the new
    # effective underlying sequence to check RBS conservation.
    orig_seq = seq_record.seq

    # If they are overlapping, see if they can be overlaid.
    # NOTE: It might seem like this should never fail, given the way that
    # _generate_permitted_feature_seq_variants() only produces variants
    # that preserve translation of the overlapped strand. However, it's
    # possible that the actual sequences for the two variants don't fit
    # together, thus we need to confirm that the two variants are actually
    # compatible on the sequence level (not just translation level).
    if upstream_feature.location.end > downstream_feature.location.start:
        # Check whether the parts of the sequences in the overlap are the
        # same in both sequences.
        upstream_overlapping_part = upstream_feature_forward_seq[
                downstream_feature.location.start -
                        upstream_feature.location.start:]
        downstream_overlapping_part = downstream_feature_forward_seq[:
                upstream_feature.location.end -
                        downstream_feature.location.start]
        can_be_overlaid = (
                upstream_overlapping_part == downstream_overlapping_part)
        if not can_be_overlaid:
            return False

        # Compute the effective underlying sequence for the RBS
        # conservation check.
        underlying_seq = (
                orig_seq[:upstream_feature.location.start] +
                upstream_feature_forward_seq +
                downstream_feature_forward_seq[
                        upstream_feature.location.end -
                                downstream_feature.location.start:] +
                orig_seq[downstream_feature.location.end:]
        )
        assert len(orig_seq) == len(underlying_seq)

    else:
        # Not overlapping, but we still prepare for RBS conservation check.
        underlying_seq = (
                orig_seq[:upstream_feature.location.start] +
                upstream_feature_forward_seq +
                orig_seq[upstream_feature.location.end:
                        downstream_feature.location.start] +
                downstream_feature_forward_seq +
                orig_seq[downstream_feature.location.end:]
        )
        assert len(orig_seq) == len(underlying_seq)

    # Check whether RBS expression is conserved.
    if (upstream_feature.strand == -1 and
            upstream_feature.type == 'CDS'):
        if not is_rbs_strength_conserved(
                original_rbs_strength_profile,
                upstream_feature,
                underlying_seq,
                rbs_allowed_error_order_of_magnitude):
            return False
    if (downstream_feature.strand == 1 and
            downstream_feature.type == 'CDS'):
        if not is_rbs_strength_conserved(
                original_rbs_strength_profile,
                downstream_feature,
                underlying_seq,
                rbs_allowed_error_order_of_magnitude):
            return False

    # Passed all tests.
    return True


def fix_same_dir_size_4_overlap(
        old_upstream_feature,
        old_downstream_feature,
        conflict_type,
        region_of_interest,
        seq_record,
        forbidden_codons,
        codon_usage_memex,
        original_rbs_strength_profile):
    """Attempts to fix a small overlap by moving the start of the
    downstream feature just far enough to resolve a conflict that
    can't be solved by synonymous swaps alone. If this works, it's likely
    better than doing the hard fix which results in copying a lot more
    sequence.

    The main target of this fix is overlaps of size 4 that might look
    something like this:

        TTATGA
          ATGACT

    which can be resolved to:
        TTATGA
             ATGACT

    The fix might be small enough where the RBS region is still in reach, but
    we check that expression as measured by the Salis RBS calculator is
    maintained.

    TODO: This is pretty jenky to catch the most common case. Figure out
        how to make it more robust.

    Returns an object with keys:
        * is_success: Whether fixing suceeded.
        * updated_seq_record: The new seq record with the changes reflected.
    """
    print '...Attempt resolve overlap of size 4...'

    fix_result = {
        'is_success': False,
    }

    # Must be facing same direction.
    if not old_upstream_feature.strand == old_downstream_feature.strand:
        return fix_result

    # Make a copy of the seq_record to make this easier for now. Maybe
    # figure out how to do this without copying later.
    new_seq_record = copy.deepcopy(seq_record)

    upstream_feature_id = old_upstream_feature.id
    upstream_feature = get_feature_by_id(
            new_seq_record, upstream_feature_id)
    assert upstream_feature

    downstream_feature_id = old_downstream_feature.id
    downstream_feature = get_feature_by_id(
            new_seq_record, downstream_feature_id)
    assert downstream_feature

    # If the features are not facing the same direction then just fail on
    # this step for now.
    if not upstream_feature.strand == downstream_feature.strand:
        return fix_result

    # Extract the forward features.
    upstream_feature_forward_seq = new_seq_record.seq[
            upstream_feature.location.start:upstream_feature.location.end]
    downstream_feature_forward_seq = new_seq_record.seq[
            downstream_feature.location.start:downstream_feature.location.end]

    if upstream_feature.strand == 1:
        assert upstream_feature_forward_seq[-1] == 'A'
    else:
        assert downstream_feature_forward_seq[0] == 'T' # Complement of A
        # Need to make the end of the upstream feature a T.
        updated_seq = (
                new_seq_record.seq[:upstream_feature.location.end - 1] +
                'T' +
                new_seq_record.seq[upstream_feature.location.end:]
        )
        assert len(new_seq_record.seq) == len(updated_seq), "%d != %d" % (
                len(new_seq_record.seq), len(updated_seq))
        new_seq_record.seq = updated_seq

    # Now copy the 3 that are overlapped, not including the last one.
    insert_seq = str(upstream_feature_forward_seq[-3:])
    insert_start_position = upstream_feature.location.end
    safe_features = [upstream_feature, downstream_feature]

    insert_sequence(
            new_seq_record,
            insert_seq,
            insert_start_position,
            safe_features=safe_features,
            insert_feature_type=None,
            insert_feature_id=None,
            insert_feature_strand=1)

    # Shift the downstream feature.
    downstream_feature = downstream_feature._shift(3)
    replace_feature(
            new_seq_record, downstream_feature.id, downstream_feature)

    # Now check that we can fix any remaining conflicts with synonymous swaps.
    synonymous_swap_result = (
        fix_conflict_using_synonymous_codon_resolution(
                upstream_feature,
                downstream_feature,
                conflict_type,
                region_of_interest,
                new_seq_record,
                forbidden_codons,
                codon_usage_memex,
                original_rbs_strength_profile
        ))
    if synonymous_swap_result['is_success']:
        # Add feature annotation for this resolution.
        feature_type = FIX_OVERLAP_SIZE_4
        feature_id = ('%s_%s_%s' % (
                feature_type,
                get_feature_gene(upstream_feature),
                get_feature_gene(downstream_feature)))
        feature_location = FeatureLocation(
                region_of_interest[0], region_of_interest[1])
        feature = SeqFeature(
                type=feature_type,
                location=feature_location,
                strand=1,
                id=feature_id
        )
        add_feature_to_seq_record(new_seq_record, feature)

        fix_result['is_success'] = True
        fix_result['updated_seq_record'] = new_seq_record
        return fix_result

    # This is not a valid fix.
    return fix_result


def fix_overlap_pair(
        upstream_feature_id, downstream_feature_id, seq_record):
    """Fixes a pair of overlapping features.

    This is the default "brute-force" fix that will always pull apart
    two genes.

    The type of fix depends on the relative strand polarity. In some cases,
    the genes may simply be pulled apart. In other cases, we additionally
    copy RBS regions for one or both overlapping genes.

    Args:
        upstream_feature_id: Id of the upstream feature.
        downstream_feature_id: Id of the downstream feature.
        seq_record: SeqRecord object that will be mutated.

    Returns:
        A Boolean indicating whether the fix was successful.
    """
    upstream_feature = get_feature_by_id(
            seq_record, upstream_feature_id)
    assert upstream_feature

    downstream_feature = get_feature_by_id(
            seq_record, downstream_feature_id)
    assert downstream_feature

    # Figure out overlap size
    overlap_size = max(0,
            upstream_feature.location.end - downstream_feature.location.start)

    print 'BEFORE HEAD COPY'
    d_seq = downstream_feature.extract(seq_record.seq)
    print d_seq[:20], '...', d_seq[-20:]


    # First split them apart.
    #     xxxxxxxxx
    #         xxxxxxxxxx
    # Becomes:
    #     xxxxxxxxx
    #              xxxxxxxxxx
    pair_overlaps = does_upstream_feature_overlap_downstream_feature(
            upstream_feature, downstream_feature)
    if pair_overlaps and overlap_size > 0:
        seq_record = shift_seq_feature_right_with_head_copy(
                downstream_feature, seq_record, upstream_feature.location.end)

    # Get the feature handles again in case they have changed. They may have
    # changed if the implementation of the split method destroyed the old
    # features and created new ones to represent moving them.
    upstream_feature = get_feature_by_id(
            seq_record, upstream_feature_id)
    downstream_feature = get_feature_by_id(
            seq_record, downstream_feature_id)

    print 'AFTER HEAD COPY'
    d_seq = downstream_feature.extract(seq_record.seq)
    print d_seq[:20], '...', d_seq[-20:]

    # Now, depending on the relative polarities, we need to do additional
    # manipulation.
    if upstream_feature.strand == 1 and downstream_feature.strand == -1:
        # Opposing case:
        #     >>>>>>>>>>
        #               <<<<<<<<<<<<
        # They are already split. Nothing more to do.
        #     >>>>>>>>>>
        #               <<<<<<<<<<<<
        pass

    elif upstream_feature.strand == 1 and downstream_feature.strand == 1:
        # Same dir,forward case:
        #     >>>Rd...Rd>>>>>>>
        #                      >>>>>>>>>>>>>
        # We need to copy the RBS site before the right one to make:
        #     >>>Rd...Rd>>>>>>>
        #                      >>Rd...Rd>>>>>>>>>>>>>>
        rbs_buffer_seq_start = (upstream_feature.location.end - overlap_size -
                RBS_BUFFER_SIZE)
        rbs_buffer_seq_end = rbs_buffer_seq_start + RBS_BUFFER_SIZE
        rbs_buffer_seq = seq_record.seq[
                rbs_buffer_seq_start:rbs_buffer_seq_end]

        insert_feature_id = (
                str(downstream_feature.id) +
                '_downstream_' +
                InsertType.FIX_OVERLAP_RBS_COPY)
        seq_record = insert_sequence(
                seq_record,
                rbs_buffer_seq,
                upstream_feature.location.end,
                safe_features=[],
                insert_feature_type=InsertType.FIX_OVERLAP_RBS_COPY,
                insert_feature_id=insert_feature_id,
                insert_feature_strand=1,
        )

    elif upstream_feature.strand == -1 and downstream_feature.strand == -1:
        # Same dir,reverse case:
        #     <<<<<<<<
        #             <<<<<Ru...Ru<<<<
        # We need to copy the RBS site before the head of the left one:
        #     <<<<<<<<<<Ru...Ru<<
        #                        <<<<<<Ru...Ru<<<<
        rbs_buffer_seq_start = downstream_feature.location.start + overlap_size
        rbs_buffer_seq_end = rbs_buffer_seq_start + RBS_BUFFER_SIZE
        rbs_buffer_seq = seq_record.seq[
                rbs_buffer_seq_start:rbs_buffer_seq_end]

        insert_feature_id = (
                str(upstream_feature.id) +
                '_upstream_' +
                InsertType.FIX_OVERLAP_RBS_COPY)
        seq_record = insert_sequence(
                seq_record,
                rbs_buffer_seq,
                upstream_feature.location.end,
                safe_features=[],
                insert_feature_type=InsertType.FIX_OVERLAP_RBS_COPY,
                insert_feature_id=insert_feature_id,
                insert_feature_strand=-1,
        )

    else:
        # Opposite dir case:
        #     <<<<Rd...Rd<<<<<<<<<
        #                         >>>>>>>>>>Ru...Ru>>>>>>
        # Need to copy both RBS sites
        # First copy downstream
        #     <<<<Ru...Ru<<<<<<<<<
        #                          >>>Rd...Rd>>>>>>>>>>>>>Rd...Rd>>>>>>
        # Then upstream
        #     <<<<Rd...Rd<<<<<<<<<<Ru..Ru<<<
        #                                   >>>Rd...Rd>> >>>>>>>>>>Ru...Ru>>>>>>

        downstream_rbs_buffer_seq_start = (upstream_feature.location.end -
                overlap_size - RBS_BUFFER_SIZE)
        downstream_rbs_buffer_seq_end = (downstream_rbs_buffer_seq_start +
                RBS_BUFFER_SIZE)
        downstream_rbs_buffer_seq = seq_record.seq[
                downstream_rbs_buffer_seq_start:downstream_rbs_buffer_seq_end]

        insert_feature_id = (
                str(downstream_feature.id) +
                '_downstream_' +
                InsertType.FIX_OVERLAP_RBS_COPY)
        seq_record = insert_sequence(
                seq_record,
                downstream_rbs_buffer_seq,
                upstream_feature.location.end,
                safe_features=[],
                insert_feature_type=InsertType.FIX_OVERLAP_RBS_COPY,
                insert_feature_id=insert_feature_id,
                insert_feature_strand=1,
        )

        # Get the right feature again since it gets moved.
        downstream_feature = get_feature_by_id(
            seq_record, downstream_feature_id)

        upstream_rbs_buffer_seq_start = (downstream_feature.location.start +
                overlap_size)
        upstream_rbs_buffer_seq_end = (upstream_rbs_buffer_seq_start +
                RBS_BUFFER_SIZE)
        upstream_rbs_buffer_seq = seq_record.seq[
                upstream_rbs_buffer_seq_start:upstream_rbs_buffer_seq_end]

        # NOTE: Also inserted at the end of the left piece so that the left
        # feature, since we're always pushing everything to the right.
        insert_feature_id = (
                str(upstream_feature.id) +
                '_upstream_' +
                InsertType.FIX_OVERLAP_RBS_COPY)
        seq_record = insert_sequence(
                seq_record,
                upstream_rbs_buffer_seq,
                upstream_feature.location.end,
                safe_features=[],
                insert_feature_type=InsertType.FIX_OVERLAP_RBS_COPY,
                insert_feature_id=insert_feature_id,
                insert_feature_strand=-1,
        )

    # Successful fix.
    return True


def _generate_permitted_feature_seq_variants(
        feature,
        other_feature,
        conflict_type,
        seq_record,
        region,
        forbidden_codons,
        codon_usage_memex,
        swap_forbidden_only=False):
    """Generates all of the sequence variants of the feature, modifying
    codons within the provided region, and not adding forbidden codons, all
    while respecting the feature that is overlapped, if relevant.
    """
    # Figure out the original sequence, accounting for polarity.
    original_feature_seq = str(feature.extract(seq_record.seq))

    # Other sequence used only when overlap.
    if conflict_type == 'overlap':
        other_feature_seq = str(other_feature.extract(seq_record.seq))
        other_feature_translation = translate_custom(other_feature_seq)

    # If the region doesn't capture the feature, just return
    # the current sequence, and no other variants.
    if not is_feature_in_region(feature, region):
        return [original_feature_seq]

    codon_indeces_to_change = get_region_codon_indeces_in_feature(
            feature, region)
    if not codon_indeces_to_change:
        return [original_feature_seq]

    def _modify_sequence(codon_index_ptr, feature_seq):
        """Recursive helper function that finds all synonymous variations
        of feature_seq.
        """
        variants = []

        codon_index = codon_indeces_to_change[codon_index_ptr]
        codon = original_feature_seq[codon_index * 3:codon_index * 3 + 3]
        assert codon, 'Bug in _modify_sequence indexing at index %d' % codon_index

        if swap_forbidden_only and not codon in forbidden_codons:
            variants.append(feature_seq)
        else:
            # Get synonymous codons, with special handling if it's the first
            # codon.
            if codon_index == 0:
                synonymous_codons = codon_usage_memex.get_start_codon_list()
            else:
                synonymous_codons = codon_usage_memex.get_synonymous_codons(codon)
            assert codon in synonymous_codons, "%s at index %d not in %s" % (
                    codon, codon_index, str(synonymous_codons))

            # Make sure the current codon is preferentially at the front
            # of the list, so that if a solution with it there works first,
            # we preferentially use it.
            synonymous_codons.remove(codon)
            synonymous_codons = [codon] + synonymous_codons

            # Loop through the syn_codons, creating new feature variants.
            for syn_codon in synonymous_codons:
                if not syn_codon in forbidden_codons:
                    new_seq = (
                            feature_seq[:codon_index * 3] +
                            syn_codon +
                            feature_seq[codon_index * 3 + 3:]
                    )
                    variants.append(new_seq)

        # If this feature overlaps the other, make sure the variant doesn't
        # change the translation of the other one.
        # This should substantially prune the tree of combinations
        # that we need to explore.
        if conflict_type == 'overlap':
            variants = filter(
                    lambda variant:
                            does_changed_feature_respect_overlapped_feature(
                                    feature,
                                    variant,
                                    other_feature,
                                    other_feature_seq,
                                    other_feature_translation
                            ),
                            variants
            )

        # If we've replaced the last codon, return the present variants.
        if codon_index_ptr == len(codon_indeces_to_change) - 1:
            return variants

        # Otherwise, make a recursive call to extend the next position.
        child_variants = []
        for variant in variants:
            child_variants.append(_modify_sequence(codon_index_ptr + 1, variant))

        # Return the flattened list of variants.
        return [child for child_list in child_variants for child in child_list]

    return _modify_sequence(0, original_feature_seq)


def _get_feature_forward_seq(feature, feature_seq):
    """Returns the forward seq of the given variant for the given
    feature.
    """
    # Get the forward strand version of both.
    if feature.strand == 1:
        return feature_seq
    else:
        return reverse_complement(feature_seq)


def does_changed_feature_respect_overlapped_feature(
        changed_feature,
        changed_feature_seq,
        overlapped_feature,
        overlapped_feature_seq,
        overlapped_feature_translation):
    """Determines whether the changing_feature_seq does not break the
    translation of the overlapped feature.
    """
    is_changed_upstream = does_upstream_feature_overlap_downstream_feature(
            changed_feature, overlapped_feature)
    if not is_changed_upstream:
        assert does_upstream_feature_overlap_downstream_feature(
                overlapped_feature, changed_feature), (
                        "Features were expected to overlap.")

    # Get the forward strand version of both.
    changed_feature_forward_seq = _get_feature_forward_seq(
            changed_feature, changed_feature_seq)
    overlapped_feature_forward_seq = _get_feature_forward_seq(
            overlapped_feature, overlapped_feature_seq)

    # Now get the new overlapped forward seq.
    # This involves nasty array math to use global positions to get local
    # sequence indeces.
    if is_changed_upstream:
        # cccccc|ccccc|
        #       |ooooo|ooooo
        overlap_size = (changed_feature.location.end -
                overlapped_feature.location.start)
        changed_seq_overlapped_part = changed_feature_forward_seq[
                -1 * overlap_size:]
        overlapped_seq_remainder = overlapped_feature_forward_seq[
                overlap_size:]
        new_overlapped_feature_forward_seq = (
                changed_seq_overlapped_part + overlapped_seq_remainder
        )
    else:
        #         |cccc|ccccccc
        # oooooooo|oooo|
        overlap_size = (overlapped_feature.location.end -
                changed_feature.location.start)
        overlapped_seq_initial = overlapped_feature_forward_seq[:
                -1 * overlap_size]
        changed_seq_overlapped_part = changed_feature_forward_seq[:
                overlap_size]
        new_overlapped_feature_forward_seq = (
                overlapped_seq_initial + changed_seq_overlapped_part
        )
    assert len(overlapped_feature_forward_seq) == len(
            new_overlapped_feature_forward_seq), "Length not conserved."

    # Correct the polarity of the new overlapped seq if necessary.
    # This puts the feature seq back into its correct frame.
    if overlapped_feature.strand == -1:
        new_overlapped_feature_seq = reverse_complement(
                new_overlapped_feature_forward_seq)
    else:
        new_overlapped_feature_seq = new_overlapped_feature_forward_seq

    # Check whether the translation is preserved.
    try:
        valid_combo = (overlapped_feature_translation ==
                translate_custom(new_overlapped_feature_seq)
        )
    except TranslationError:
        # One situation where this error is thrown is when the start codon
        # is not a valid start codon in the table.
        valid_combo = False

    return valid_combo


# DEBUGGING

if __name__ == '__main__':
    from biopython_util import get_genome_record
    from refactor_config import ORIGINAL_GENOME_RECORD
    genome_record = get_genome_record('../data/mds42_full.gbk')
    with open(FIX_OVERLAPS_CACHE_LOCATION, 'w') as fh:
        pickle.dump(ORIGINAL_GENOME_RECORD, fh)
