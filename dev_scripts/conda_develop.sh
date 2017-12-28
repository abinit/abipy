#!/bin/bash
set -ev  # exit on first error, print each command

echo Installing AbiPy dependencies with conda.
echo Adding conda-forge, matsci and abinit to channels
echo Working in CONDA_PREFIX: ${CONDA_PREFIX}
conda config --add channels conda-forge
conda config --add channels matsci
conda config --add channels abinit

# https://github.com/ContinuumIO/anaconda-issues/issues/542
conda install -y -c anaconda setuptools
conda install -y --file ./requirements.txt
conda install -y --file ./requirements-optional.txt
# Install bader (http://theory.cm.utexas.edu/henkelman/code/bader/) from matsci
conda install -y -c matsci bader

# Install abinit from abinit conda channel.
conda install -y -c gmatteo abinit=8.6.1
abinit --version
abinit --build

echo Installation complete. Use: conda install abinit to install Fortran executable