"""
Script for extracting a region of a genome and all its features and inserting
them into another genome.
"""

import copy

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


def excise_genome_region(genbank_record_or_path, region_interval,
        output_path=None, feature_transform_fn_list=[],
        output_format='genbank'):
    """Extract all features contained within a region of a genome
    and write the result to a new output file.

    Only features that are completely contained by the region will
    be included.

    Args:
        genbank_record_or_path: File location of the source genome,
            or a SeqRecord.
        region_interval: Two-tuple giving the 0-indexed region to excise.
        output_path: If provided, the extracted region will be written to this
            output file.
        feature_transform_fn_list: List of functions to apply to every
            feature in the resulting genome.

    Returns:
        A SeqRecord representation of the excised region with positions zeroed
        relative to the interval.
    """
    if isinstance(genbank_record_or_path, SeqRecord):
        genome_record = copy.deepcopy(genbank_record_or_path)
    else:
        genome_record = SeqIO.read(genbank_record_or_path, 'genbank')

    # Find and extrac the features contained in the region.
    extracted_feature_list = []
    for feature in genome_record.features:
        if _is_feature_contained_by_interval(feature, region_interval):
            extracted_feature_list.append(feature)

    # Adjust the feature positions relative to the interval.
    transformed_feature_list = [feature._shift(-1 * region_interval[0])
            for feature in extracted_feature_list]

    # Apply transform functions to features.
    for transform_fn in feature_transform_fn_list:
        transformed_feature_list = [transform_fn(feature) for feature in
                transformed_feature_list]

    # Build the new seq record to contain the modified features.
    excised_seq_record = SeqRecord(genome_record.seq[region_interval[0]:
            region_interval[1]])
    excised_seq_record.features = transformed_feature_list

    # Write the result.
    if output_path:
        SeqIO.write(excised_seq_record, output_path, output_format)

    return excised_seq_record


def _is_feature_contained_by_interval(feature, interval):
    """Checks whether a feature is completely inside the interval.

    Args:
        feature: SeqFeature object.
        interval: Two-tuple representing pythonic interval.

    Returns:
        A boolean.
    """
    return (feature.location.start >= interval[0] and
            feature.location.end <= interval[1])


def overlay_features_in_genome_at_position(feature_source,
        target_genome, target_start_position, output_path):
    """
    Args:
        feature_source: Source path or SeqRecord containing the features
            to overlay on the target.
        target_genome: Target genome path or SeqRecord.
        target_start_position: Pythonic start position in the target genome.
        output_path: Path to write the resulting Genbank file to.
    """
    if isinstance(feature_source, SeqRecord):
        source_genome_record = feature_source
    else:
        raise ValueError("Invalid type for feature_source.")

    target_genome_record = SeqIO.read(target_genome, 'genbank')

    # Shift the source features relative to the target start position.
    transformed_feature_list = [feature._shift(target_start_position)
            for feature in source_genome_record.features]

    # Overlay the features on the target.
    target_genome_record.features.extend(transformed_feature_list)

    SeqIO.write(target_genome_record, output_path, 'genbank')


if __name__ == '__main__':
    def remove_unwanted_qualifiers(feature):
        UNWANTED_QUALIFIERS = ['modified_by', 'created_by', 'ref', 'alt']
        for qualifier in UNWANTED_QUALIFIERS:
            if qualifier in feature.qualifiers:
                del feature.qualifiers[qualifier]
        return feature

    lambda_prophage_seq_record = excise_genome_region(
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/C321(F1,M1,Ao)_I3_not_strict.gb',
            (805241, 817465),
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/lambda_prophage.genbank',
            [remove_unwanted_qualifiers]
            )

    overlay_features_in_genome_at_position(lambda_prophage_seq_record,
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/rec1_c321d.genbank',
            805241,
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/rec1_c321d_with_lambda_prophage_annotation.genbank')
