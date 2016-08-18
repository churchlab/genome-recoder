"""
Global context object that stores handles to reusable objects
and data sources.

NOTE: Previously we were using the refactor_config module for this purpose
but it has grown unwieldy at this point. So for now we have a bit of a messy
hybrid.
"""

import os
import pickle

from codon_replacer import GraphSearchCodonReplacer
from feature_profile import CodonRarityFeatureProfile
from feature_profile import GCContentFeatureProfile
from feature_profile import SecondaryStructureFeatureProfile
from refactor_config import CACHE_DIR
from refactor_config import CODONS_TO_REMOVE
from refactor_config import ORIGINAL_CODON_USAGE_MEMEX
from refactor_config import REFACTORED_CODON_USAGE_MEMEX


PROFILE_VALUES_CACHE = os.path.join(
        CACHE_DIR, 'mds42_profile_values.cache')


class RefactorContext(object):
    """Object that stores handles to reusable objects and data sources.

    NOTE: This object is not guaranteed to be thread-safe.
    """

    def __init__(self, genome_record, use_cache_data=True):
        """Constructor.
        """
        self.set_genome_record(genome_record)

        self.use_cache_data = use_cache_data

        self._ss_feature_profile_factory = (
                SecondaryStructureFeatureProfile.factory())

        self._forbidden_codon_set = set(CODONS_TO_REMOVE)

        # The following are lazily initiliazed the first time they are called.
        # This makes sure that any setters can be more safely used by the
        # client as these objects are often interdependent.
        self._codon_replacer_obj = None
        self._feature_id_to_profile_values_map = None


    def get_genome_record(self):
        """Returns the SeqRecord representing the genome.

        Clients should be careful about mutating this.
        """
        return self._genome_record


    def set_genome_record(self, genome_record):
        """Returns the SeqRecord representing the genome.

        Clients should be careful about mutating this.
        """
        self._genome_record = genome_record


    def get_ss_feature_profile_factory(self):
        """Returns the SecondaryStructureFeatureProfile factory.

        NOTE: This object should not be shared among different threads.
        """
        return self._ss_feature_profile_factory


    def get_codon_replacer_obj(self):
        """Returns the codon replacer object.
        """
        if not self._codon_replacer_obj:
            self._codon_replacer_obj = GraphSearchCodonReplacer(
                    self.get_forbidden_codon_list(),
                    self.get_original_codon_usage_memex())
        return self._codon_replacer_obj


    def get_feature_id_to_profile_values_map(self, force_recalculate=False):
        """Returns a map from feature_id to a map keyed by profile
        name of the original values array for each.

        This map gets created lazily.
        """
        if self._feature_id_to_profile_values_map is None:
            self._feature_id_to_profile_values_map = (
                    calc_original_feature_profile_values(
                        self._genome_record,
                        self._ss_feature_profile_factory,
                        force_recalculate=force_recalculate,
                        use_cache_data=self.use_cache_data))
        return self._feature_id_to_profile_values_map


    def set_feature_id_to_profile_values_map(self, feature_id_to_profile_values_map):
        """Explicitly set the map.
        """
        self._feature_id_to_profile_values_map = (
                feature_id_to_profile_values_map)


    def get_forbidden_codon_set(self):
        """Returns the list of codons we are removing in this refactor.
        """
        return self._forbidden_codon_set


    def set_forbidden_codon_set(self, forbidden_codon_set):
        """Set the forbidden codon set for this object as per a builder pattern.
        """
        assert isinstance(forbidden_codon_set, set)
        self._forbidden_codon_set = forbidden_codon_set


    def get_forbidden_codon_list(self):
        """Returns the list of codons we are removing in this refactor.
        """
        return list(self.get_forbidden_codon_set())


    def get_original_codon_usage_memex(self):
        """Returns the codon usage object with the original usage profile.
        """
        return ORIGINAL_CODON_USAGE_MEMEX


    def get_refactored_codon_usage_memex(self):
        """Returns the codon usage object with the original usage profile.
        """
        return REFACTORED_CODON_USAGE_MEMEX



###############################################################################
# Utility methods
###############################################################################

def calc_original_feature_profile_values(
        genome_record, ss_feature_profile_factory, force_recalculate=False,
        use_cache_data=True):
    """Returns a map from feature_id to a map keyed by profile
    name of the original values array for each.
    """
    print '...Getting original feature profiles...'

    if use_cache_data:
        if not os.path.exists(CACHE_DIR):
            os.mkdir(CACHE_DIR)

        if not force_recalculate:
            try:
                with open(PROFILE_VALUES_CACHE) as fh:
                    feature_id_to_profile_values_map = pickle.load(fh)
                    print ('......Using cached profile values from %s.' %
                            PROFILE_VALUES_CACHE)
                    return feature_id_to_profile_values_map
            except IOError:
                # No cache.
                pass

    feature_id_to_profile_values_map = {}

    features_to_profile = filter(
            lambda feature: feature.type == 'CDS',
            genome_record.features
    )

    codon_rarity_profile_kwargs = {
            'original_codon_usage_memex': ORIGINAL_CODON_USAGE_MEMEX,
            'refactored_codon_usage_memex': REFACTORED_CODON_USAGE_MEMEX
    }

    num_features = len(features_to_profile)
    for feature_index in range(num_features):
        feature = features_to_profile[feature_index]
        feature_id = feature.id

        # Debug output.
        print 'Calculating original profile for feature %d of num_features %d' % (
                    feature_index + 1, num_features)
        print '...Feature id: %s' % feature_id

        # Construct the FeatureProfile objects.
        gc_content_profile = GCContentFeatureProfile(
                feature_id, genome_record)

        secondary_structure_profile = ss_feature_profile_factory(
                feature_id, genome_record)

        codon_rarity_profile = CodonRarityFeatureProfile(
                feature_id, genome_record, **codon_rarity_profile_kwargs)

        feature_id_to_profile_values_map[feature_id] = {
                gc_content_profile.__class__.get_name():
                        gc_content_profile.values,
                secondary_structure_profile.__class__.get_name():
                        secondary_structure_profile.values,
                codon_rarity_profile.__class__.get_name():
                        codon_rarity_profile.values,
        }

    if use_cache_data:
        # Cache the results.
        with open(PROFILE_VALUES_CACHE, 'w') as fh:
            pickle.dump(feature_id_to_profile_values_map, fh)

    return feature_id_to_profile_values_map


if __name__ == '__main__':
    from biopython_util import get_genome_record
    genome_record = get_genome_record('../data/mds42_full.gbk')
    rc = RefactorContext(genome_record)
    rc.get_feature_id_to_profile_values_map(True)
