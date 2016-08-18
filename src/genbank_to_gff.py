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
