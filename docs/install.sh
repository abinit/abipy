#!/bin/bash
set -ev  # exit on first error, print each command

pip install -r ./requirements.txt
pip install -r ../requirements-optional.txt
pip install -r ../requirements-panel.txt
pip install -r requirements.txt
conda install graphviz -c conda-forge --yes
