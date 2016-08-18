#!/usr/bin/env python

"""Script to add copyrights to our source files.
"""

import os


COPYRIGHT_TEXT = (
    "# Copyright (C) 2016 The President and Fellows of Harvard College\n"
    "#\n"
    "# This file may require additional software or modules to be installed to run\n"
    "# run properly.\n"
    "#\n"
    "# The Genome Recoder software is available under an internal non-commercial\n"
    "# research and academic use license.  Questions about this software or the\n"
    "# licensing thereof can be addressed to  Office of Technology Development,\n"
    "# Harvard University, email: otd@harvard.edu.\n"
    "#\n"
    "# @author Gleb Kuznetsov (kuznetsov@g.harvard.edu)\n"
    "\n"
)

for dirpath, dirnames, filenames in os.walk('src'):
    if 'rbs_calc' in dirpath:
        continue
    for f in filenames:
        if os.path.splitext(f)[1] == '.py':
            file_path = os.path.join(dirpath, f)
            with open(file_path) as fh:
                contents = fh.read()
            contents = COPYRIGHT_TEXT + contents
            with open(file_path, 'w') as output_fh:
                output_fh.write(contents)
