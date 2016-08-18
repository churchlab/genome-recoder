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
Script to convert .genbank file to GFF3 format.
"""

from Bio import SeqIO
from BCBio import GFF

def convert(genbank_in_path, gff_out_path):
    """Converts a genbank file to a GFF file.
    """
    with open(genbank_in_path) as in_handle:
        with open(gff_out_path, 'w') as out_handle:
            GFF.write(SeqIO.parse(in_handle, 'genbank'), out_handle)


if __name__ == '__main__':
    convert(
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/rec1_c321d.preliminary.gbk',
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/rec1_c321d.preliminary.gbk.gff'
    )
