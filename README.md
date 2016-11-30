# rE.coli-57 software

Software library for recoding genomes used in **Ostrov et al. (2016). Design, synthesis, and testing toward a 57-codon genome. Science, 353(6301)**.

Please contact us if you're interested in using this software for your own recoding projects: Gleb Kuznetsov (kuznetsov@g.harvard.edu)

## Setup

### Setup python environment and install python dependencies

The software was developed and tested using python 2.7.

The required python libraries are listed in `requirements.txt`. We recommend
using [virtualenv](http://pypi.python.org/pypi/virtualenv) to create a
sandboxed environment.

    virtualenv venv

    source venv/bin/activate .

    pip install -r requirements.txt

### Install Unafold (required for secondary structure calculation).

1. Download and install from here:

    <http://homepages.rpi.edu/~zukerm/download/UNAFold_download.html>

If using Linux, you probably want to download the RPM and use these instructions:
<https://overlappingminds.com/sh/thoughts/0aa6d79e-fb8b-4f84-b287-f8e4494eac49>

### Install Salis Ribosome Binding Site (RBS) Calculator

    cd src/

    git clone git@github.com:hsalis/Ribosome-Binding-Site-Calculator-v1.0.git rbs_calc

    cd rbs_calc/

    touch __init__.py

### Install NUPACK (Nucleic Acid Package) - Required by Salis Calculator

1. Download here [http://www.nupack.org/](http://www.nupack.org/). Our software
   expects Nupack to be located at /opt/nupack/.  So either install it there, or update `src/refactor_config.py`.

2. Copy the files in bin/ to /usr/local/bin

3. Test that `hybrid-ss-min` can be run from the shell. Test the following command in a new terminal:

        hybrid-ss-min -h

## Running tests

If using virtualenv, be sure the correct python environment is activated.

    cd src/

    ./run_tests.sh

## Usage

The main entry point to the code is `src/main.py`. After setting up the proper Python environment, you should be able to simply run:

    cd src/

    python main.py

Results are written to `src/output`.

## License

Copyright (C) 2016 The President and Fellows of Harvard College

The Genome Recoder software is available under an internal non-commercial
research and academic use license.  Questions about this software or the
licensing thereof can be addressed to  Office of Technology Development,
Harvard University, email: otd@harvard.edu.

See the LICENSE file.

