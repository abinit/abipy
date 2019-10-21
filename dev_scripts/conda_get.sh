#!/bin/bash
set -e  # exit on first error

# Install conda with travis: https://conda.io/docs/travis.html
if [[ "${TRAVIS_OS_NAME}" == "osx" ]]; then
    curl -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh;
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b -p ${HOME}/miniconda
export PATH="${HOME}/miniconda/bin:${PATH}"

hash -r

conda config --set always_yes yes --set changeps1 yes
conda update -q conda
# Useful for debugging any issues with conda
conda info -a
conda config --add channels conda-forge
echo "Installing abinit from abinit channel in abinit-environment..."
conda create -q -n abinit-environment python=${TRAVIS_PYTHON_VERSION}
source activate abinit-environment
conda install -y -c abinit abinit=${ABINIT_VERSION}
abinit --version
abinit --build
echo "Creating test-environment for python stack..."
conda create -q -n test-environment python=${TRAVIS_PYTHON_VERSION}
source activate test-environment
