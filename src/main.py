"""
Main entry point for performing genome refactoring.

NOTE(2/6/13): Started encapsulating strategy logic in strategies.py.
"""

from refactor_config import DEBUG
from refactor_config import ORIGINAL_GENOME_RECORD
from strategies import refactor_with_min_overlap_fixes_and_preserve_rbs


def main(tmp_file_prefix):
    """Entry point into performing refactoring.

    NOTE: There is not really any command line argument support so for now
    the primary way for clients to configure parameters is to do edit code in
    refactor_config.
    """
    print 'Refactoring genome ...'

    genome_record = ORIGINAL_GENOME_RECORD

    # Update the name of the genome record so that it's easier to distinguish
    # from the original in Geneious.
    genome_record.name += '_mod'

    # Run the refactoring strategy.
    refactor_with_min_overlap_fixes_and_preserve_rbs(
            genome_record,
            tmp_file_prefix,
            debug=DEBUG)


if __name__ == '__main__':
    import cProfile
    from datetime import datetime
    import os

    from refactor_config import OUTPUT_DIR

    TMP_FILE_PREFIX = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    CPROFILE_OUTPUT_DEST = os.path.join(
            OUTPUT_DIR, TMP_FILE_PREFIX + '_cprofile.out')
    cProfile.run('main(TMP_FILE_PREFIX)', CPROFILE_OUTPUT_DEST)
