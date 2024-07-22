#!/bin/bash
set -ev  # exit on first error, print each command

pip install -r requirements.txt
conda install graphvix -c conda-forge --yes
