#!/bin/band
# WARNING: Must be executed on a linux box.
conda env remove -n abipy_binder
conda create -n abipy_binder python=3.8
python -m pip install --editable "../[optional]"
conda install abinit -c conda-forge
conda env export > environment.yml
