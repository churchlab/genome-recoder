"""
Object that maintains a score profile across a feature, per codon. This is used
to constrain choices of synonymous codons while refactoring. Further, choosing
synonymous codons may be constrained by multiple profiles at the same time.

NOTE: We are in the process of figuring out exactly how this will work.
"""

import os
import pickle
import subprocess

from Bio.SeqUtils import GC

import biopython_util
from biopython_util import get_genome_record
from tool_wrappers.hybrid_ss_min_wrapper import HybridSSMinWrapper


class FeatureProfile(object):
    """Base class for a Profile. Subclasses should override with their own
    logic for computing scores, etc.

    We store a datapoint in the scores member variable, one for each codon.
    This datapoint is calculated depending on the rule for type of
    FeatureProfile (subclasses specify their own).  Generally this will
    be some function of a sliding window ending on that codon. This strategy
    then lends itself to a graph-search across synonymous codon substitutions.
    """

    # The feature this profile is for.
    original_seq_record = None

    # The feature extracted from original_seq_record.
    # Cache the lookup/convenience reference.
    feature = None

    # The original (polarity-aware) sequence for the feature
    original_feature_seq = None

    # Underlying sequence for the entire genome, accounting for feature polarity
    # so that calculations are done in the correct relative direction.
    polarity_aware_underlying_seq = None

    # Start location of the feature, taking into account polarity.
    polarity_aware_feature_location_start = None

    # End location of the feature, taking into account polarity.
    polarity_aware_feature_location_end = None

    # The list of values, indexed by codon index in the feature.
    # (i.e. values[0] is te profile value for the first codon in the feature.
    values = None

    # Tuple representing the bounds for the values.
    value_range = (0, 1.0)

    # How much error is allowed between a profile and a proposed alternative.
    error_tolerance = 0.0

    # Amount to increment the error_tolerance by in case we need looser bounds
    # to find a valid solution.
    error_inc = 0.0

    # Size of the window in the sliding window strategy that most
    # FeatureProfiles use.
    sliding_window_size = None

    # Stores results of computations in case they are repeated.
    computation_cache = None

    failed_error_treshold_count = None

    def __init__(
            self,
            feature_id,
            original_seq_record,
            **kwargs):
        """Default constructor.

        Args:
            feature_id: The string id of the feature to profile.
            original_seq_record: The SeqRecord object with that contains all the
                data.
            kwargs: Dictionary of any of the following keys
                * error_tolerance
                * error_inc
                * value_range
        """
        self.original_seq_record = original_seq_record
        self.feature = biopython_util.get_feature_by_id(
                self.original_seq_record, feature_id)
        self.original_feature_seq = str(
                self.feature.extract(self.original_seq_record.seq))
        self.computation_cache = {}
        self.failed_error_treshold_count = 0

        if 'error_tolerance' in kwargs:
            self.error_tolerance = kwargs['error_tolerance']

        if 'error_inc' in kwargs:
            self.error_inc = kwargs['error_inc']

        if 'value_range' in kwargs:
            self.value_range = kwargs['value_range']

        if 'sliding_window_size' in kwargs:
            self.sliding_window_size = kwargs['sliding_window_size']

        if self.feature.strand == 1:
            self.polarity_aware_underlying_seq = self.original_seq_record.seq
            self.polarity_aware_feature_location_start = (
                    self.feature.location.start)
            self.polarity_aware_feature_location_end = (
                    self.feature.location.end)
        elif self.feature.strand == -1:
            self.polarity_aware_underlying_seq = (
                    self.original_seq_record.seq.reverse_complement())
            total_len = len(self.polarity_aware_underlying_seq)
            self.polarity_aware_feature_location_start = (
                    total_len -  self.feature.location.end)
            self.polarity_aware_feature_location_end = (
                    self.polarity_aware_feature_location_start +
                    len(self.feature))
        else:
            raise ValueError("No feature strand.")

        # Some basic validation.
        self._validate()

        # Set or compute the original score profile.
        if 'values' in kwargs:
            num_codons = len(self.feature) / 3
            self.values = kwargs['values']
            assert num_codons == len(self.values), ("Wrong number of "
                    "values passed to FeatureProfile constructor.")
        else:
            self.compute_original_feature_values()


    def _validate(self):
        """Basic validation to be performed during instantiation.
        """
        # Validation to make sure sub-classes are properly extended.
        assert hasattr(self, '_name'), "Subclasses must define _name"

        # We expect a coding feature to have an integral number
        # of codons (sets of 3-bases).
        error_msg = "Expected feature length to be multiple of 3"
        assert len(self.feature) % 3 == 0, error_msg


    def compute_original_feature_values(self):
        """Computes the profile score for each codon in the feature.
        """
        num_codons = len(self.feature) / 3
        self.values = [
                self.compute_codon_value(
                        codon_index,
                        feature_seq=self.original_feature_seq,
                        computing_original_profile=True)
                for codon_index in range(num_codons)]


    def compute_codon_value(self,
            codon_index, feature_seq, computing_original_profile=False):
        """Prepares for computing, acting as a multiplexer for determining
        the sequence to compute over, before calling the internal method
        (generally overriden by child classes), which contains the computation
        logic.
        """
        return self.compute_codon_value_internal(
                codon_index, feature_seq, computing_original_profile)


    def compute_codon_value_internal(self,
            codon_index, feature_seq, computing_original_profile=False):
        """Default computation strategy.

        Sub-classes should override this.
        """
        return 1.0


    def get_codon_value_for_codon_index(self, codon_index):
        """Returns the profile value for the given codon_index.
        """
        return self.values[codon_index]


    def get_window_seq(self, codon_index, feature_seq):
        """Return the sequence window for the computation.
        """
        # The codon start position relative to the original SeqRecord.
        codon_start_pos = self.codon_index_to_seq_record_pos(codon_index)

        # Calculating starting position of the window.
        # Graphically, if the codon we are about it TTT:
        #      NNNNNNNNNNNNNNNNNNNNNNTTT
        #     |<---- window_size ------>|
        window_start = codon_start_pos - self.sliding_window_size + 3
        window_start = max(0, window_start)

        # Figure out the starting position within the feature.
        first_pos_in_feature_seq = max(0,
                window_start - self.polarity_aware_feature_location_start)

        # And get the last position within the feature
        window_end = window_start + self.sliding_window_size
        last_pos_in_feature_seq = (window_end -
                self.polarity_aware_feature_location_start)

        # Assemble the window sequence.
        window_seq = ''
        if first_pos_in_feature_seq == 0:
            window_seq += self.polarity_aware_underlying_seq[
                    window_start:self.polarity_aware_feature_location_start]
        window_seq += feature_seq[
                first_pos_in_feature_seq:last_pos_in_feature_seq]
        assert self.sliding_window_size == len(window_seq), (
                "Expected window size: %d. Observed: %d" % (
                        self.sliding_window_size, len(window_seq)))
        return str(window_seq)


    def codon_index_to_seq_record_pos(self, codon_index):
        """Returns the start position in the seq record for
        the codon index of the given feature.
        """
        return self.polarity_aware_feature_location_start + codon_index * 3


    def calc_feature_score_relative_to_profile(
            self, codon_index, feature_seq):
        """Returns a normalized error score for the feature given its sequence,
        at the codon_index position, relative to the feature profile.

        Subclasses should override this.

        Args:
            codon_index: An integer specifying which codon in the feature
                we are scoring (generally this is the last codon in a sliding
                window).
            feature_seq: The altered sequence of the feature we are evaluating.

        Returns:
            An object with keys:
                * error_score: float
                * passes_error_threshold: Boolean
        """
        original_value = self.get_codon_value_for_codon_index(codon_index)

        altered_feature_value = self.compute_codon_value(
                codon_index, feature_seq)

        error_score = altered_feature_value - original_value
        passes_error_threshold = abs(error_score) <= self.error_tolerance
        if not passes_error_threshold:
            self.failed_error_treshold_count += 1

        return {
                'profile_name': self.get_name(),
                'error_score': error_score,
                'passes_error_threshold': passes_error_threshold
        }


    def increment_error_tolerance(self):
        """Increments the error tolerance of the profile.

        This can be used, for example, when a search over the possible space
        of sequence solutions gets into a loop or is having trouble
        progressing otherwise.

        TODO: Should there be a maximum tolerance at which point we report
        an error while trying to refactor that particular genome feature?
        """
        self.error_tolerance += self.error_inc


    def reset_failed_error_threshold_count(self):
        """Resets the failed threshold count to 0."""
        self.failed_error_treshold_count = 0


    @classmethod
    def get_name(cls):
        """Returns the name of the class.

        One place this is used is as a key into the dictionary of
        FeatureProfile values.
        """
        return cls._name



class GCContentFeatureProfile(FeatureProfile):
    """FeatureProfile that profiles the GC content across a gene.
    """

    _name = 'gc_content_profile'

    # The size of the window over which the GC content is calculated.
    sliding_window_size = 40

    error_tolerance = 0.1

    error_inc = 0.1


    def compute_codon_value_internal(self,
            codon_index, feature_seq, computing_original_profile=False):
        """Computes the GC content of the sequence with the codon
        as the last element in the sequence, and including the previous
        sliding_window_size number of bases.

        Args:
            codon_index: The index of the codon in the feature.
            source_seq: The feature sequence to computer over.

        TODO: We can make this more efficient by leveraging calculation
        from previous frame rather than recalculating the entire frame.
        """
        window_seq = self.get_window_seq(codon_index, feature_seq)
        if window_seq in self.computation_cache:
            return self.computation_cache[window_seq]
        gc_content = GC(window_seq) / 100
        self.computation_cache[window_seq] = gc_content
        return gc_content



class SecondaryStructureFeatureProfile(FeatureProfile):
    """FeatureProfile that profiles the secondary structure across
    a genome feature.
    """

    _name = 'ss_profile'

    # The size of the window over which secondary structure is calculated.
    sliding_window_size = 40

    error_tolerance = 5.0

    error_inc = 5.0

    # Object that wraps the logic for computing secondary structure.
    ss_computer = None

    def __init__(
            self,
            feature_id,
            original_seq_record,
            **kwargs):
        """Override the default constructor in order to instantiate the
        object that computes the secondary structure free energy score.
        """
        super(SecondaryStructureFeatureProfile, self).__init__(
                feature_id, original_seq_record, **kwargs)

        if 'ss_computer' in kwargs:
            self.ss_computer = kwargs['ss_computer']
        else:
            self.ss_computer = HybridSSMinWrapper()


    def compute_codon_value_internal(self,
            codon_index, feature_seq, computing_original_profile=False):
        """Computes the free energy for the window ending at the
        given codon index of the source_seq.

        Args:
            codon_index: The index of the codon in the feature.
            source_seq: The feature sequence to computer over.
        """
        window_seq = self.get_window_seq(codon_index, feature_seq)
        if window_seq in self.computation_cache:
            return self.computation_cache[window_seq]
        delta_g = self.ss_computer.compute_delta_g(window_seq)
        self.computation_cache[window_seq] = delta_g
        return delta_g


    def compute_original_feature_values(self):
        """Manually override to test if passing all the sequences to
        hybrid-ss-min is possible, and faster.
        """
        num_codons = len(self.feature) / 3

        window_seqs = [
                self.get_window_seq(codon_index, self.original_feature_seq)
                for codon_index
                in range(num_codons)
        ]

        cmd = "hybrid-ss-min --NA=RNA -q " + ' '.join(window_seqs)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        all_results = p.stdout.read()
        separated = all_results.split('\n')
        self.values = map(
                lambda x: float(x.strip()),
                filter(lambda y: bool(y), all_results.split('\n')))

        # Cache the computations.
        for window_seq, value in zip(window_seqs, self.values):
            self.computation_cache[window_seq] = value


    @classmethod
    def factory(cls):
        """Factory method that allows reuse of the HybridSSMinWrapper
        if desired.
        """
        ss_computer = HybridSSMinWrapper()
        def create_instance_fn(feature_id, original_seq_record, **kwargs):
            kwargs['ss_computer'] = ss_computer
            return cls(feature_id, original_seq_record, **kwargs)
        return create_instance_fn


class CodonRarityFeatureProfile(FeatureProfile):
    """FeatureProfile that preserves codon rarity.
    """

    _name = 'codon_rarity_profile'

    error_tolerance = 0.2

    error_inc = 0.1

    # 5 codons
    sliding_window_size = 15

    # Assigned in the constructor.
    original_codon_usage_memex = None
    refactored_codon_usage_memex = None

    def __init__(self, feature_id, original_seq_record, **kwargs):
        """Override the default constructor to require a CodonUsageMemex.

        This constructor expects a codon_usage_memex in **kwargs.
        """
        assert 'original_codon_usage_memex' in kwargs
        assert 'refactored_codon_usage_memex' in kwargs

        self.original_codon_usage_memex = kwargs['original_codon_usage_memex']
        self.refactored_codon_usage_memex = kwargs['refactored_codon_usage_memex']

        # Call the super constructor last since the initial profiling
        # requires the codon_usage_memex.
        super(CodonRarityFeatureProfile, self).__init__(
                feature_id, original_seq_record, **kwargs)


    def compute_codon_value_internal(self,
            codon_index, feature_seq, computing_original_profile=False):
        """Computes the average rarity of the codons in the window.
        """
        window_seq_str = self.get_window_seq(codon_index, feature_seq)
        if window_seq_str in self.computation_cache:
            return self.computation_cache[window_seq_str]

        # Split the window into a list of codons.
        codons = [window_seq_str[start:start + 3]
                        for start in range(0,len(window_seq_str), 3)]

        # Choose the correct codon_usage_memex
        if computing_original_profile:
            codon_usage_memex = self.original_codon_usage_memex
        else:
            codon_usage_memex = self.refactored_codon_usage_memex

        # Get usage of each codon in the window.
        codon_usage = map(
                lambda codon: self.original_codon_usage_memex.get_codon_usage(
                        codon),
                codons)

        # Return the average usage across the window.
        avg_codon_usage = sum(codon_usage) / len(codon_usage)
        self.computation_cache[window_seq_str] = avg_codon_usage
        return avg_codon_usage


class WeissmanFeatureProfile(FeatureProfile):
    """FeatureProfile that requires that works to preserve the Weissman
    ribosomal pausing profile.

    According to Weissman, ribosomal pausing is caused by RBS-like sequences.
    Previous hypotheses have argued that codon rarity accounts for ribosomal
    pausing.

    The way the current implementation will work is that anywhere that there is
    an RBS "peak" as defined by propensity to anneal with the
    anti-Shine-Dalgarno sequence (??CCTCCT) as defined in Weissman's paper,
    we want to preferentially make sure that forbidden codons are replaced so a
    to maximize annealing at that position with this test sequence.

    We believe that this profile's vote should generally be stronger than other
    profiles' votes at places where a Shine-Dalgarno-like sequence is detected.
    Otherwise, it should not report any error.
    """

    _name = 'weissman_profile'

    def compute_codon_value_internal(
            self, codon_index, feature_seq, computing_original_profile=False):
        # TODO: Implement.
        return 0
