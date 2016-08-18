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
Module for storing constants for paths relevant to anlaysis.

Previously this was inside refactor_config.py but, it makes sense to put these
constants in a separate module for at least 2 reasons:
    * Modules often need data paths but not the other heavy loading
        refactor_config does at start up.
    * We are generally trying to break up refactor_config.py
"""

import os


PWD = os.path.dirname(os.path.realpath(__file__ ))

DATA_DIR = os.path.join(PWD, '../data')

CONFIG_DIR = os.path.join(PWD, '../config')

EC_ALIGNMENT_DATA_DIR = os.path.join(PWD, '../ec_alignments')

GENOMES_DIR = os.path.join(DATA_DIR, 'genomes')

OUTPUT_DIR = os.path.join(PWD, 'output')
if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)

TMP_DATA_DIR = os.path.join(PWD, 'tmp')
if not os.path.exists(TMP_DATA_DIR):
    os.mkdir(TMP_DATA_DIR)

CACHE_DIR = os.path.join(PWD, 'cache')
if not os.path.exists(CACHE_DIR):
    os.mkdir(CACHE_DIR)

MG1655_SOURCE = os.path.join(GENOMES_DIR, 'mg1655', 'mg1655.genbank')
