#!/bin/band
# WARNING: Must be executed on a linux box.
conda env remove -n abipy_binder
conda create -n abipy_binder python=3.8
conda install --file ../requirements.txt -y -c conda-forge
conda install --file ../requirements-optional.txt -y -c conda-forge
conda install abinit -c conda-forge
conda env export > environment.yml
