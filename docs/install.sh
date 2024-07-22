#!/bin/bash
set -ev  # exit on first error, print each command

pip install -r requirements.txt
conda install -y --file ./requirements.txt
conda install graphviz -c conda-forge --yes
pip install -r ../requirements-optional.txt
pip install -r ../requirements-panel.txt
