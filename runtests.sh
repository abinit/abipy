#!/bin/bash
set -ev  # exit on first error, print each command

abicheck.py
nosetests abipy -v --with-coverage --cover-package=abipy --logging-level=INFO
# This is to run the integration tests (slow)
#; pytest -v abipy/integration_tests

# Documentation
./docs/install.sh
cd ./docs && make && cd ../

