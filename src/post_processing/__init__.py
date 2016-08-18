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
Package that contains modules for performing processing after the
main codon replacement. These might be additional constraints imposed by
synthesis requirements.

This might include actoins such as:
    * removing restriction enzyme sites
    * removing homopolymer runs.
    * normalizing regions of extreme GC.
"""
