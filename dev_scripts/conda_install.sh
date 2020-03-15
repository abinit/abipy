#!/bin/bash
set -e  # exit on first error

echo "Installing bader executable (http://theory.cm.utexas.edu/henkelman/code/bader/) from matsci ..."
conda install -y -c matsci bader
conda install -c abinit apscheduler==2.1.0
#pip install apscheduler==2.1.0

echo "Installing requirements listed requirements.txt and requirements-optional.txt ..."
# https://github.com/ContinuumIO/anaconda-issues/issues/542
conda install -y -c anaconda setuptools
conda install nomkl
conda install -y --file ./requirements.txt
conda install -y --file ./requirements-optional.txt
conda install -y -c conda-forge graphviz python-graphviz

echo "Installation completed"
