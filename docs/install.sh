#!/bin/bash
set -ev  # exit on first error, print each command

python -m pip install --editable "../[optional,panel,docs]"
conda install graphviz -c conda-forge --yes
