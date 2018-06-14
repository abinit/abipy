#!/bin/bash
set -e  # exit on first error

# Install conda with travis: https://conda.io/docs/travis.html
if [[ "${TRAVIS_OS_NAME}" == "osx" ]]; then
    if [[ "${TRAVIS_PYTHON_VERSION}" == "2.7" ]]; then
        curl -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh;
    else
        curl -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh;
    fi
else
    if [[ "${TRAVIS_PYTHON_VERSION}" == "2.7" ]]; then
        wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh;
    else
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    fi
fi

bash miniconda.sh -b -p ${HOME}/miniconda
export PATH="${HOME}/miniconda/bin:${PATH}"

hash -r

conda config --set always_yes yes --set changeps1 yes
conda update -q conda
# Useful for debugging any issues with conda
conda info -a
# Replace dep1 dep2 ... with your dependencies
conda create -q -n test-environment python=${TRAVIS_PYTHON_VERSION}
source activate test-environment

