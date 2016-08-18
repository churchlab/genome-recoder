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
Utility mehods for replacing codons.

The reason for this module is that we ended up with circular dependencies
due to requiring both RefactorContext and CodonReplacer objects. The quickest
way to resolve this for now is to create this new module.
"""

from datetime import datetime
import os
import pickle
import multiprocessing
import stat

from biopython_util import update_seq_record_feature
from feature_profile import CodonRarityFeatureProfile
from feature_profile import GCContentFeatureProfile
from feature_profile import SecondaryStructureFeatureProfile
from refactor_context import RefactorContext


def get_essential_feature_ids(genome_record):
    """Returns the list of feature ids that we want to change codons for.
    """
    FEATURES_TO_REFACTOR = ['CDS']
    return [feature.id for feature in genome_record.features
            if feature.type in FEATURES_TO_REFACTOR]


def replace_forbidden_codons(
        refactor_context,
        num_cores=1,
        tmp_file_prefix=None,
        debug=False):
    """Distribute the task of fixing features over the number of cores
    provided.

    This can actually be completely parallelized since we are no longer
    changing feature start/end locations, but only the sequences themselves.
    It can be done safely if features that are near each other are done
    sequentially, so that they take into account changes to the other
    when adhering to a profile. However, it's probably safe to run them
    completely independent of each other, as the feature profiles will
    be preserved between even any adjacent features, even if not exact.
    """
    print 'Replacing codons ...'
    genome_record = refactor_context.get_genome_record()
    codons_to_remove = refactor_context.get_forbidden_codon_list()
    original_codon_usage_memex = (
            refactor_context.get_original_codon_usage_memex())
    refactored_codon_usage_memex = (
            refactor_context.get_refactored_codon_usage_memex())

    TMP_PARTIAL_DATA_DIR = 'tmp'

    # Make sure the directory for partial data exists.
    if not os.path.exists(TMP_PARTIAL_DATA_DIR):
        os.mkdir(TMP_PARTIAL_DATA_DIR)
        # User and group have all permissions | get group id from directory.
        os.chmod(TMP_PARTIAL_DATA_DIR, stat.S_ISGID | 0775)

    # Get the features that we'll fix.
    essential_feature_ids = get_essential_feature_ids(genome_record)

    # Partition the input depending on the number cores.
    partition_size = len(essential_feature_ids) / num_cores
    partition_start_positions = [i * partition_size for i in range(num_cores)]
    partition_end_positions = [start + partition_size for start in
            partition_start_positions]
    # Manually fix the last end position to end exactly at the number of
    # features.
    partition_end_positions[-1] = len(essential_feature_ids)

    if not tmp_file_prefix:
        tmp_file_prefix = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')

    feature_id_to_profile_values_map = (
            refactor_context.get_feature_id_to_profile_values_map())

    # Track processes
    processes = {}

    for core in range(num_cores):

        range_start = partition_start_positions[core]
        range_end = partition_end_positions[core]
        p_id = '%s_partial_%d_%d' % (tmp_file_prefix, range_start, range_end)
        tmp_result_file = os.path.join(
                TMP_PARTIAL_DATA_DIR, p_id + '.pickle')

        # TODO: For performance, investigate whether we should be passing
        # copies of the objects passed as args to the child process. None
        # of these are mutated internally, so it's safe as we are doing now,
        # but perhaps shared memory might be wonky.
        p = multiprocessing.Process(
                target=replace_codons_in_feature_subset,
                args=(
                        genome_record,
                        essential_feature_ids,
                        codons_to_remove,
                        original_codon_usage_memex,
                        refactored_codon_usage_memex,
                        feature_id_to_profile_values_map,
                        range_start,
                        range_end,
                        tmp_result_file,
                        debug
                ))
        p.start()
        processes[p_id] = {
                'process': p,
                'partial_results_file': tmp_result_file,
        }


    for p_id, process_obj in processes.items():
        process_obj['process'].join()

    # We'll also return some metadata in case we want to do anything with it.
    metadata = {}

    # Replace all the codon features in the genome_record
    # once all the subsets are computed.
    print "Combining results and updating genome record..."
    for p_id, process_obj in processes.items():
        update_seq_record_with_partial_results(
                genome_record, metadata, process_obj['partial_results_file'])
    print "...Done."

    # Return the genome_record with all the features updated.
    return genome_record, metadata


def replace_codons_in_feature_subset(
        genome_record,
        essential_feature_ids,
        codons_to_remove,
        original_codon_usage_memex,
        refactored_codon_usage_memex,
        feature_id_to_profile_values_map,
        range_start,
        range_end,
        tmp_result_file,
        debug):
    """Work on a subset of the features to fix.
    """
    assert range_end >= range_start

    results = {}

    refactor_context = RefactorContext(genome_record)
    refactor_context.set_feature_id_to_profile_values_map(
        feature_id_to_profile_values_map)

    # Code path that helps debug across a limited range when debug is True.
    effective_range_start = range_start
    effective_range_end = range_end if not debug else range_start + 1

    num_features = len(essential_feature_ids)
    for feature_index in range(effective_range_start, effective_range_end):
        print 'Fixing feature %d of  %d' % (feature_index + 1, num_features)
        feature_id = essential_feature_ids[feature_index]
        print 'Feature id: %s' % feature_id

        result = replace_codons_in_single_feature(refactor_context, feature_id)

        # Add the results to the growing dictionary.
        results[feature_id] = result

    # Write resuts to file as soon as we're done.
    with open(tmp_result_file, 'w') as fh:
        pickle.dump(results, fh)


def replace_codons_in_single_feature(
        refactor_context,
        feature_id,
        explicit_genome_record=None,
        start_codon_index=0,
        last_codon_index=None,
        avoid_codons_in_positions={}):
    """Method that encapsulates creating all the objects necessary to
    perform codon replacement and then calls the appropriate methods of
    those objects to complete replacement.

    The underlying genome_record in refactor_context is not mutated.
    """
    # Pull what we need from refactor_context.
    if explicit_genome_record:
        genome_record = explicit_genome_record
    else:
        genome_record = refactor_context.get_genome_record()
    codon_replacer_obj = refactor_context.get_codon_replacer_obj()
    feature_id_to_profile_values_map = (
            refactor_context.get_feature_id_to_profile_values_map())
    profile_values = feature_id_to_profile_values_map[feature_id]
    ss_feature_profile_factory = (
            refactor_context.get_ss_feature_profile_factory())
    original_codon_usage_memex = (
            refactor_context.get_original_codon_usage_memex())
    refactored_codon_usage_memex = (
            refactor_context.get_refactored_codon_usage_memex())

    # Construct the FeatureProfile objects.
    gc_content_profile_kwargs = {
            'values': profile_values[GCContentFeatureProfile.get_name()]
    }
    gc_content_profile = GCContentFeatureProfile(
            feature_id, genome_record, **gc_content_profile_kwargs)
    # print '...Generated gc_content_profile'

    secondary_structure_profile_kwargs = {
        'values': profile_values[SecondaryStructureFeatureProfile.get_name()]
    }
    secondary_structure_profile = ss_feature_profile_factory(
            feature_id, genome_record,
            **secondary_structure_profile_kwargs)
    # print '...Generated secondary_structure_profile'

    codon_rarity_profile_kwargs = {
            'original_codon_usage_memex': original_codon_usage_memex,
            'refactored_codon_usage_memex': refactored_codon_usage_memex,
            'values': profile_values[CodonRarityFeatureProfile.get_name()],
    }
    codon_rarity_profile = CodonRarityFeatureProfile(
            feature_id, genome_record, **codon_rarity_profile_kwargs)
    # print '...Generated codon_rarity_profile'

    # Modify the profile list used by the CodonReplacer.
    feature_profile_list = [
            gc_content_profile,
            secondary_structure_profile,
            codon_rarity_profile
    ]
    codon_replacer_obj.set_profiler_list(feature_profile_list)

    # Run the replacement algorithm.
    result = codon_replacer_obj.replace_codons_in_feature(
            feature_id,
            genome_record,
            start_codon_index=start_codon_index,
            last_codon_index=last_codon_index,
            avoid_codons_in_positions=avoid_codons_in_positions)

    return result


def update_seq_record_with_partial_results(
        seq_record, metadata, pickled_partial_results_loc):
    """Update the SeqRecord object with the pickled results.

    Args:
        seq_record: SeqRecord to be mutated.
        metadata: Dictionary of metadata to be updated.
        pickled_partial_results_loc: Relative location of the pickled
            partial result object.
    """
    try:
        with open(pickled_partial_results_loc) as fh:
            results = pickle.load(fh)
            metadata.update(results)
            for feature_id, result_obj in results.iteritems():
                if result_obj['is_success']:
                    update_seq_record_feature(
                            seq_record,
                            feature_id,
                            result_obj,
                    )
    except IOError:
        # Results file may not exists.
        pass
