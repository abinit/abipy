#!/bin/bash
set -e  # exit on first error

echo "Installing AbiPy dependencies with conda."
echo "Adding conda-forge and abinit to channels"
echo "Working in CONDA_PREFIX: ${CONDA_PREFIX} ..."
conda config --add channels conda-forge

echo "Installing abinit from abinit channel before pymatgen..."
conda install -y -c abinit abinit=${ABINIT_VERSION}
abinit --version
abinit --build

echo "Installing bader executable (http://theory.cm.utexas.edu/henkelman/code/bader/) from matsci ..."
conda install -y -c matsci bader apscheduler==2.1.0
pip install apscheduler==2.1.0

echo "Installing requirements listed requirements.txt and requirements-optional.txt ..."
# https://github.com/ContinuumIO/anaconda-issues/issues/542
conda install -y -c anaconda setuptools
#conda install nomkl
conda install -y --file ./requirements.txt
conda install -y --file ./requirements-optional.txt
conda install -y -c conda-forge graphviz python-graphviz

echo "Testing abinit after pymatgen installation ..."
#echo "Installing abinit from abinit channel ..."
#conda install -y -c abinit abinit=${ABINIT_VERSION}
abinit --version
abinit --build

echo "Installation completed"
