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
Methods for interfacing with regulonDB data.
"""

import os

from numpy import array
import pandas as pd


# Directory containing this script
PWD = os.path.dirname(os.path.realpath(__file__ ))

# Data location
DATA_DIR = os.path.join(PWD, '../data')

# Location of regulonDB-related data.
REGULON_DB_DIR = os.path.join(DATA_DIR, 'regulondb')

# Operons dataset.
REGULON_DB_OPERONS = os.path.join(REGULON_DB_DIR, 'OperonSet_no_comments.txt')

# Genes dataset. We specifically use the gene products dataset where we don't
# need the sequences.
REGULON_DB_GENE_PRODUCTS = os.path.join(REGULON_DB_DIR,
        'GeneProductSet_no_comments.txt')

# Keio knockout study data.
KEIO_ESSENTIAL_ONLY = os.path.join(REGULON_DB_DIR, 'keio_clean_essential_only.csv')

# Destination for output report file.
OUTPUT_REPORT_FILE = os.path.join(REGULON_DB_DIR, 'operon_essentiality_report.csv')


def analyze_operon_essential_genes():
    """Joins operon data to essential gene data.

    Returns:
        The DataFrame containing operons data.
    """
    # We start out by carrying out a basic analysis in a procedural fashion.
    # We can break out separate functions as we go. The basic task is to join
    # three datasets: operons, genes, and essentiality, and output a table
    # reporting the aspects of interest in this join, including:
    #     * number of essential genes per operon
    #     * size of the operon
    #     * operon location
    #     * essential genes in the operon

    # First read in relevant data.
    operons_df = pd.DataFrame.from_csv(REGULON_DB_OPERONS, sep='\t')
    keio_essential_df = pd.DataFrame.from_csv(KEIO_ESSENTIAL_ONLY)

    # Get the set of all essential genes.
    keio_essential_gene_set = set(keio_essential_df.gene)

    # Add column related to essential genes.
    def _compute_essential_genes_column(comma_separated_genes):
        gene_list = comma_separated_genes.split(',')
        essential_gene_list = [gene for gene in gene_list if gene in
                keio_essential_gene_set]
        return ','.join(essential_gene_list)
    operons_df['essential_genes'] = operons_df.genes.map(
            _compute_essential_genes_column)
    operons_df['num_essential_genes'] = operons_df.essential_genes.map(
            lambda x: len(x) and len(x.split(',')))

    # Other useful columns.
    operons_df['operon_size'] = (operons_df.end.map(lambda x: int(x)) -
            operons_df.start.map(lambda x: int(x)) + 1)

    # Write the result
    operons_df.to_csv(OUTPUT_REPORT_FILE, sep='\t')

    return operons_df


def identify_mg1655_partition_positions(partition_size_lower_bound=40000,
        partition_size_upper_bound=52000):
    """Identify regions in MG1655 that are safe to partition at. That is, these
    are regions that lie outside of operons.

    Args:
        partition_size_lower_bound: Lower bound for the size of partitions we
            are interested.
        partition_size_upper_bound: Upper bound for the same.
    """
    # Read in the RegulonDB data.
    operons_df = pd.DataFrame.from_csv(REGULON_DB_OPERONS, sep='\t')
    genes_df = pd.DataFrame.from_csv(REGULON_DB_GENE_PRODUCTS, sep='\t')
    genes_df = genes_df.set_index('gene_name')

    # Sort the dataframe by starting point.
    operons_df = operons_df.sort(columns='start')

    # Add a column that takes the difference between consecutive and sometimes
    # overlapping intervals.  I'm a pandas noob so there is probably a better
    # way to do the delta than what follows.
    end_points = array([0] + list(operons_df['end']))
    start_points = array(list(operons_df['start']) + [0])
    delta_previous_interval = start_points - end_points
    operons_df['delta_previous_interval'] = delta_previous_interval[:-1]

    # Now we want to extract in-between intervals. These will be intervals
    # in between non-overlapping operons. The bounds of the interval are
    # inclusive.
    # To allow for getting the corresponding genes below, we want
    # inbetween_interval list to have one-to-one correspondence with
    # the dataframe, so we store None when there is no in-between interval
    # that preceeds the current operon.
    inbetween_intervals = [None]
    for i in range(1, len(operons_df)):
        delta_previous_interval = operons_df.iloc[i][
                'delta_previous_interval']
        if delta_previous_interval > 0:
            start = operons_df.iloc[i]['start']
            interval = (start - delta_previous_interval + 1, start - 1)
            inbetween_intervals.append(interval)
        else:
            inbetween_intervals.append(None)
    assert len(inbetween_intervals) == len(operons_df)

    # Now that we have in-between intervals, we want to find those that
    # surround regions of the size that we are interested in with our
    # partitions.

    # As a first step, for each in-between interval, we find:
    #     * set of in-between intervals that follow it after a jump that falls
    #         within the partition size.
    #     * the gene immediately following the interval.
    interval_map = {}
    num_inbetween_intervals = len(inbetween_intervals)
    for current_idx, current_interval in enumerate(inbetween_intervals):
        if current_interval is None:
            continue

        interval_map[current_interval] = {}

        # Identify the first gene in the operon immediately following this
        # interval.
        gene_list = operons_df.iloc[current_idx]['genes'].split(',')

        # Create list of tuples of genes with their positions.
        genes_with_positions = []
        for gene in gene_list:
            regulondb_gene = genes_df.loc[gene]
            genes_with_positions.append((
                    regulondb_gene.name,
                    int(regulondb_gene['start']),
                    int(regulondb_gene['end'])))
            # HACK: (as if there weren't any others) Add the bnumber as a gene.
            genes_with_positions.append((
                    regulondb_gene['bnumber'],
                    int(regulondb_gene['start']),
                    int(regulondb_gene['end'])))

        # Sort the genes by start position.
        genes_with_positions = sorted(genes_with_positions,
                key=lambda x: x[1])

        # Add the genes to the map.
        interval_map[current_interval]['genes'] = genes_with_positions

        # Find intervals following this one.
        interval_map[current_interval]['following_intervals'] = []
        next_interval_idx = current_idx + 1
        while next_interval_idx < num_inbetween_intervals:
            next_interval = inbetween_intervals[next_interval_idx]
            next_interval_idx += 1
            if next_interval is None:
                continue
            gap = next_interval[0] - current_interval[1]
            if gap < partition_size_lower_bound:
                # Not big enough gap, check next.
                continue
            if gap > partition_size_upper_bound:
                # All following are too far away.
                break
            # If passing the above gauntlet, then this next interval is in the
            # valid jump range.
            interval_map[current_interval]['following_intervals'].append(
                    next_interval)

    return interval_map


if __name__ == '__main__':
    identify_mg1655_partition_positions()
