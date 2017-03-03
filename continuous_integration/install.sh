#!/usr/bin/env bash
# This script is inspired from **pgmpy** implementation of continous test
# integration. This is meant to "install" all the packages required for installing
# semantic.

# License: The MIT License (MIT)

set -e

# sudo apt-get update -qq
# sudo apt-get install build-essential -qq

if [[ "$DISTRIB" == "conda" ]]; then
	# Deactivate the travis-provided virtual environment and setup a
	# conda-based environment instead
	deactivate

	# Use the miniconda installer for faster download / install of conda
	# itself
	wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh \
		-O miniconda.sh
    chmod +x miniconda.sh && ./miniconda.sh -b
    export PATH=~/miniconda2/bin:$PATH
	conda config --set always_yes yes --set changeps1 no
	conda update conda
	conda info -a

	conda create -n testenv python=$PYTHON_VERSION --file continuous_integration/requirements.txt
    source activate testenv
fi

if [[ "$COVERAGE" == "true" ]]; then
	pip install coverage coveralls
fi

# Build pytftb
python setup.py develop
