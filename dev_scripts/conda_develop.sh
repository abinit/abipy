#!/bin/bash
set -ev  # exit on first error, print each command

echo Installing abipy dependencies with conda.
echo Adding conda-forge, mastci and abinit to channels
echo Working in CONDA_PREFIX: %{CONDA_PREFIX}
conda config --add channels conda-forge
conda config --add channels matsci
conda config --add channels abinit

conda install -y --file ./requirements.txt
conda install -y --file ./requirements-optional.txt

echo Installation complete. Use: conda install abinit to install Fortran executable