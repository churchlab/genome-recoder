"""
Functions for generating vendor orders.

VOCAB:
    * segment: Delivery unit (50kb as of 7/25/14)
    * fragment: Synthesis unit (3kb for Gen9 as of 7/25/14)
"""

import csv
import re

import pandas as pd

from collections import defaultdict
from biopython_util import get_feature_label
from biopython_util import InsertType
from post_processing.gen9_synthesis_constraints_checker import check_all
from well_id_generator import WellIdGenerator


GEN9_ORDER_CSV_FIELD_NAMES = [
    'Plate',
    'Well',
    'Customer_ID',
    'Sequence',
]


def create_gen9_order(
        genome_record,
        order_csv_path,
        start_seg_num=None,
        end_seg_num=None,
        upstream_homology_arm=None, downstream_homology_arm=None,
        upstream_flanking_cut_site=None, downstream_flanking_cut_site=None,
        force_seg_start_in_first_col_of_plate=False,
        start_plate=1,
        start_well_letter='A',
        start_well_number=1):
    """Gen9 accepts a csv with a row per [3kb] fragment.

    Args:
        genome_record: SeqRecord.
        order_csv_path: Output csv path.
        start_seg_num: Lowest segment number to add.
        end_seg_num: Highest segment number to add.
        upstream_homology_arm: If provided, the upstream arm for the whole
            [50kb] segment (not each fragment). Forward strand.
        downstream_homology_arm: If provided, the downstream arm for the whole
            [50kb] segment (not each fragment). Forward strand.
        upstream_flanking_cut_site: If provided, added to every fragment.
            Must also provide downstream_flanking_cut_site.
        downstream_flanking_cut_site: Downstream analog of upstream.
        force_seg_start_in_first_col_of_plate: If True, start each segment in
            first column of plate. Write blank rows elsewhere.
    """
    if upstream_flanking_cut_site:
        assert downstream_flanking_cut_site

    well_id_generator = WellIdGenerator(
            start_well_letter=start_well_letter,
            start_well_number=start_well_number,
            start_plate=start_plate,
            include_plate=True)

    frag_features = [feature for feature in genome_record.features
            if feature.type == InsertType.PARTITION_GEN9_SEG]

    # For the first and last fragment, we want to add homology to pYES1L.
    # We separate the fragment features by segment id so we can identify
    # the first and last fragment in each segment.
    seg_id_to_frag_feature_map = defaultdict(list)
    for frag_feature in frag_features:
        frag_label = get_feature_label(frag_feature)
        seg_id = re.match(r'seg[0-9]+', frag_label).group()
        seg_num = int(re.search(r'[0-9]+', seg_id).group())
        if ((start_seg_num is not None and seg_num < start_seg_num) or
                (end_seg_num is not None and seg_num > end_seg_num)):
            continue
        seg_id_to_frag_feature_map[seg_id].append(frag_feature)
    seg_id_to_frag_feature_map__sorted = {}
    for seg_id in seg_id_to_frag_feature_map.iterkeys():
        seg_id_to_frag_feature_map__sorted[seg_id] = sorted(
            seg_id_to_frag_feature_map[seg_id],
            key=get_feature_label)

    with open(order_csv_path, 'w') as csv_fh:
        writer = csv.DictWriter(csv_fh, GEN9_ORDER_CSV_FIELD_NAMES)
        writer.writeheader()
        for seg_id in sorted(seg_id_to_frag_feature_map.iterkeys()):
            sorted_frag_feature_list = seg_id_to_frag_feature_map__sorted[
                    seg_id]

            added_upstream_homology_arm = False
            added_downstream_homology_arm = False
            for idx, frag_feature in enumerate(sorted_frag_feature_list):
                plate, well_id = well_id_generator.next()

                # If forcing segments to align with columns, skip
                # rows if necessary.
                if idx == 0 and force_seg_start_in_first_col_of_plate:
                    while not re.search('01', well_id):
                        writer.writerow({
                            'Plate': plate,
                            'Well': well_id,
                            'Customer_ID': '',
                            'Sequence': '',
                        })
                        plate, well_id = well_id_generator.next()

                fragment_id = get_feature_label(frag_feature)

                fragment_seq = frag_feature.extract(genome_record.seq)

                # Maybe add homology arms.
                if upstream_homology_arm is not None and idx == 0:
                    fragment_seq = upstream_homology_arm + fragment_seq
                    added_upstream_homology_arm = True
                if (downstream_homology_arm is not None and
                        idx == len(sorted_frag_feature_list) - 1):
                    fragment_seq += downstream_homology_arm
                    added_downstream_homology_arm = True

                # Maybe add flanking cut sites.
                if upstream_flanking_cut_site:
                    assert downstream_flanking_cut_site
                    fragment_seq = (
                            upstream_flanking_cut_site +
                            fragment_seq +
                            downstream_flanking_cut_site)

                writer.writerow({
                    'Plate': plate,
                    'Well': well_id,
                    'Customer_ID': fragment_id,
                    'Sequence': str(fragment_seq),
                })

            # Make sure we added the pYES1L flanking pieces.
            if upstream_homology_arm:
                assert added_upstream_homology_arm
            if downstream_homology_arm:
                assert added_downstream_homology_arm


def verify_gen9_order(order_csv, verbose=False):
    """Verify that the order doesn't violate synthesis constraints.

    Returns:
        Dict mapping from seg_id to verification result.
    """
    result = {}
    gen9_fragments_df = pd.read_csv(order_csv)
    for idx, row in gen9_fragments_df.iterrows():
        seg_id = row['Customer_ID']
        if pd.isnull(seg_id):
            continue
        if verbose:
            print '\n\n>>>>>>>>>>>', seg_id
        seq = row['Sequence']
        seq_verification_result = check_all(seq)
        result[seg_id] = seq_verification_result
        if verbose:
            for key in seq_verification_result.iterkeys():
                if len(seq_verification_result[key]):
                    print seq_verification_result[key]
    return result
