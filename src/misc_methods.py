"""
Catch all for manual methods that don't quite fit in any other module.

Some of the stuff might be old but we're saving it just in case
it becomes useful.
"""


def update_genome_record_from_existing_partials(tmp_file_prefix):
    """NEEDS FIXING

    Special script developed during debugging that lets us skip
    unnecessarily re-calculating gene refactorings, by taking the
    partial output files containing pickled refactored genomes and directly
    updating the genome record after fixing overlaps (also probaly from cache).
    """
    GENOME_SOURCE = '../data/mds42_full.gbk'

    # Get the genome record.
    genome_record = get_genome_record(GENOME_SOURCE)

    # Update the name of the genome record so that it makes sense in display.
    genome_record.name += '_mod'

    aa_to_codons_dict = build_codon_usage_dict()
    original_codon_usage_memex = CodonUsageMemex(aa_to_codons_dict)

    # Fix overlaps.
    # Fix overlaps.
    cfp = ConflictingPairFixer(
            genome_record,
            CODONS_TO_REMOVE,
            original_codon_usage_memex,
            cache=True,
    )
    genome_record = cfp.fix_overlaps()

    # Data from refactoring that will be pickled, and can then be explored
    # pythonically for analysis or debugging purposes.
    metadata = {}

    partial_files = filter(
            lambda path: re.match(tmp_file_prefix, path),
            os.listdir('tmp'))

    for partial in partial_files:
        partial_path = os.path.join('tmp', partial)
        update_seq_record_with_partial_results(
                genome_record, metadata, partial_path)

    # Write the resulting genome.
    genome_output_file = os.path.join(
            OUTPUT_DIR, tmp_file_prefix + '_mds42_refactored.gbk')
    with open(genome_output_file, 'w') as output_fh:
        SeqIO.write(genome_record, output_fh, 'genbank')

    # Write the metadata.
    metadata_output_file = os.path.join(
            OUTPUT_DIR, tmp_file_prefix + '_mds42_refactored.meta')
    with open(metadata_output_file, 'w') as metadata_output_fh:
        pickle.dump(metadata, metadata_output_fh)
