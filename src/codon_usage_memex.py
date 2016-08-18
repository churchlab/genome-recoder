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
Utility for working with codon usage distributions.
"""

import csv
import os
import random

PWD = os.path.dirname(os.path.realpath(__file__ ))

STOP_CODON_AMINO_ACID_SYMBOL = '*'

# Codons with usage below this threshold are considered rare.
RARITY_THRESHOLD = 0.15

class CodonUsageMemex:
    """Object that wraps a codon usage table, providing
    utility methods such as getting synonymous codons
    for a codon, getting amino acids, etc.
    """

    aa_to_codons_dict = None
    codon_to_aa_dict = None

    # E-coli specific. May want to check this if doing other species.
    # NOTE: Ordered by rarity.
    start_codons = ['ATG', 'GTG', 'TTG', 'ATT']

    stop_codons = None


    def __init__(self, aa_to_codons_dict):
        """Constructor.

        NOTE: The factory constructor classmethod from_source_file() will
        usually be a better way to create an object of this type.
        """
        self.aa_to_codons_dict = aa_to_codons_dict
        self._build_codon_to_aa_dict()

        # Cache the stop codons since we use these as a sanity check
        # when reading through codons.
        if STOP_CODON_AMINO_ACID_SYMBOL in self.aa_to_codons_dict:
            self.stop_codons = set(
                    self.aa_to_codons_dict.get(
                            STOP_CODON_AMINO_ACID_SYMBOL, {}).keys())
        else:
            self.stop_codons = set()


    def _build_codon_to_aa_dict(self):
        """Make reverse dictionary, codon to amino acid.
        """
        self.codon_to_aa_dict = {}
        for aa, codon_dict in self.aa_to_codons_dict.iteritems():
            for codon_name in codon_dict.iterkeys():
                self.codon_to_aa_dict[codon_name] = aa


    def get_synonymous_codons(self, codon, include_forbidden=True,
            this_codon_first=False):
        """Returns the list of synonymous.

        We want to avoid biasing codons based on order, so we return a random
        ordering each time, biased by the projected codon usage.

        NOTE: This may seriously slow things down, so look closely at this
        to make sure it still works.
        """
        amino_acid = self.codon_to_aa_dict[codon]
        ordered_codons = self.aa_to_codons_dict[amino_acid].keys()

        # The list we will return.
        randomized_order = []

        # First split between remaining and forbidden codons.
        remaining_codons = []
        forbidden_codons = []
        for this_codon in ordered_codons:
            usage = self.get_codon_usage(this_codon)
            if this_codon_first and this_codon == codon:
                randomized_order.append(this_codon)
            elif usage == 0:
                forbidden_codons.append(this_codon)
            else:
                remaining_codons.append(this_codon)

        # Sample one at a time, with bias proportional to target codon usage.
        while len(remaining_codons) > 0:
            if len(remaining_codons) == 1:
                randomized_order.append(remaining_codons[0])
                remaining_codons = []
                break

            # Get the usage sum for normalization.
            codon_usage_sum = 0
            for codon in remaining_codons:
                codon_usage_sum += self.get_codon_usage(codon)

            # Build a normalized cumulative ordering.
            cum_score = 0
            codons_with_scores = []
            for codon in remaining_codons:
                cum_score += self.get_codon_usage(codon) / codon_usage_sum
                codons_with_scores.append((codon, cum_score))
            codons_with_scores[-1] = (codons_with_scores[-1][0], 1.0)

            # Sample the next one.
            boom = random.random()
            found = False
            for codon, score in codons_with_scores:
                if boom <= score:
                    randomized_order.append(codon)
                    remaining_codons.remove(codon)
                    found = True
                    break
            assert found

        if include_forbidden:
            randomized_order += forbidden_codons

        return randomized_order


    def get_start_codon_list(self):
        return self.start_codons[:]


    def is_stop_codon(self, codon):
        """Checks whether a codon is a stop codon.
        """
        return codon in self.stop_codons


    def get_codon_usage(self, codon):
        """Returns a float between 0 and 1 indicating how often the codon
        is used in coding sequences.

        Codon usage has been shown to have effects at least on translation,
        in addition to being conserved across species, and so we take this
        measure into account when doing synonymous codon swapping.
        """
        try:
            amino_acid = self.codon_to_aa_dict[codon]
            codon_data = self.aa_to_codons_dict[amino_acid][codon]
            return codon_data['usage']
        except KeyError:
            return 0


    @classmethod
    def from_source_file(cls, source_file_location):
        """Factory method for creating a CodonUsageMemex from a source file.
        """
        # TODO: Implement.
        pass


    @classmethod
    def build_from_removed_codons_usage_report(cls,
            removed_codons_usage_report_location):
        """Factory method that builds a CodonUsageMemx using the data
        in the codon usage report created in
        mds42_analysis.analyze_codon_usage and then with a column
        called 'projected_usage' manually added to ensure rarity.

        Returns a CodonUsageMemex instance.
        """
        aa_to_codons_dict = {}
        with open(removed_codons_usage_report_location) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                # Parse the data.
                amino_acid = row['amino_acid']
                codon = row['codon']
                attempted_remove = row['attempted_remove'] == 'TRUE'
                if attempted_remove:
                    codon_usage = 0
                else:
                    codon_usage = row['projected_usage']
                    if not codon_usage:
                        codon_usage = row['original_usage']
                    codon_usage = float(codon_usage)

                # Make sure an entry exists for the amino acid.
                if not amino_acid in aa_to_codons_dict:
                    aa_to_codons_dict[amino_acid] = {}

                # Write the relevant data.
                aa_to_codons_dict[amino_acid][codon] = {
                        'usage': codon_usage
                }

        # Now go back through and make sure the usages add up to one,
        # only adjusting non-rare codons, as defined by RARITY_THRESHOLD,
        # when the distribution is off.
        for amino_acid, codons_dict in aa_to_codons_dict.iteritems():
            balance_codon_usage(amino_acid, codons_dict)

        return CodonUsageMemex(aa_to_codons_dict)


###############################################################################
# Helper Methods
###############################################################################

DEFAULT_CODON_USAGE_SOURCE_FILE = os.path.join(PWD,
        '../data/ecoli-codon-usage.txt')

def build_codon_usage_dict(
        codon_usage_source_file=DEFAULT_CODON_USAGE_SOURCE_FILE):
    """Build codon and reverse codon table with only the codons that will
    be used.

    For each amino acid, we get the list of codons that code for that
    amino acid with their associated distributions. Each codon is represented
    as a list, with the first element being

    Returns a two-tuple of dictionaries:
        * Amino acid to list of codons map for.
        * The reverse map, codon to amino acid.
    """
    # Read the source file.
    data_as_strings = []
    with open(codon_usage_source_file) as fh:
        for line in fh:
            splitline = line.split()
            for i in [0, 5, 10, 15]:
                data_as_strings.append(splitline[i:i+5])

    # Build a dictionary from amino acid to list of codons.
    aa_to_codons_dict = {}
    for string in data_as_strings:
        dna_codon = string[0].replace('U','T')
        amino_acid = string[1]
        codon_usage = float(string[2])
        if not amino_acid in aa_to_codons_dict:
            aa_to_codons_dict[amino_acid] = {}
        aa_to_codons_dict[amino_acid][dna_codon] = {
                'usage': codon_usage
        }

    return aa_to_codons_dict
    
def build_usage_dict_from_gbk(genome,
    codon_usage_source_file='../data/ecoli-codon-usage.txt'):
    """Build codon and reverse codon table from a seqrecord object,
    instead of from a static codon count table. 
    
    This allows us to get codon usage from other organisms or genome 
    variants.
    
    However, we still start with the table in order to get the 
    codon/amino acid associations.
    """
    
    codon_lookup = CodonUsageMemex(build_codon_usage_dict(
        codon_usage_source_file))
    codon_usage = {}
    aa_count = {}
    
    for codon, aa in codon_lookup.codon_to_aa_dict.items():
        codon_usage[codon] = codon_lookup.aa_to_codons_dict[aa][codon]
        codon_usage[codon]['count'] = 0
        codon_usage[codon]['usage'] = 0
        aa_count[aa] = 0
        
        
    for feature in genome.features:
        
        #skip non-gene features
        if feature.type != 'CDS': continue
        gene_start = feature.location.start.position
        gene_end = feature.location.end.position
        
        seq = genome[gene_start:gene_end].seq
        if feature.strand < 0: seq = seq.reverse_complement()
        seq = seq.tostring()
        if len(seq) % 3 is not 0: continue
        
        #split gene into codons, add up counts
        for i in range(0,len(seq)-3):
            if i % 3 != 0: continue
            codon = seq[i:i+3]
            aa = codon_lookup.codon_to_aa_dict[codon]
            codon_usage[codon]['count'] += 1
            aa_count[aa] += 1
            
    #after counting is complete, find usage
    for codon, aa in codon_lookup.codon_to_aa_dict.items():
        codon_usage[codon]['usage'] = float(codon_usage[codon]['count']) / float(aa_count[aa])
        codon_usage[codon]['RSCU'] = codon_usage[codon]['usage'] * len(
            codon_lookup.aa_to_codons_dict[aa])
    
    return(codon_lookup.aa_to_codons_dict)

def balance_codon_usage(amino_acid, codons_dict):
    """Utility methods that updates the codon usage values for the particular
    amino acid so that the total comes to 1.0.
    """
    EPSILON = 0.05

    usage_sum = sum(map(
            lambda codon_dict: codon_dict['usage'],
            codons_dict.values()))
    usage_undershoot = 1.0 - usage_sum

    # Rough accounting for float error.
    # NOTE: It's actually not necessary that we get the usage absolutely
    # correct, given the "heuristic" nature of the whole feature profile
    # strategy for refactoring.
    if ((usage_undershoot > 0 and usage_undershoot - EPSILON <= 0) or
            (usage_undershoot < 0 and usage_undershoot + EPSILON >= 0)):
        # Nothing to do. The usage error is not significant.
        return

    # First identify the non-rare codons to redistribute to.
    # NOTE: We intentionally want to keep a few codons disproportionately rare.
    codons_to_update = []
    for codon, data in codons_dict.iteritems():
        if data['usage'] > RARITY_THRESHOLD:
            codons_to_update.append(codon)

    # Now redistribute the amount proportionately.
    usage_denominator = sum([
            codons_dict[codon]['usage'] for codon in codons_to_update])
    for codon in codons_to_update:
        codon_share = float(codons_dict[codon]['usage']) / usage_denominator
        codons_dict[codon]['usage'] += usage_undershoot * codon_share


if __name__ == '__main__':
    # codon_usage_dict = build_codon_usage_dict()
    # for codon, data in codon_usage_dict.iteritems():
    #     print codon, data

    CodonUsageMemex.build_from_removed_codons_usage_report(
            '../data/mds42_codon_usage_analysis_with_projections.csv')

