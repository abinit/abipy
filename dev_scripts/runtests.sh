#!/bin/bash
set -ev  # exit on first error, print each command

abicheck.py

nosetests abipy -v --with-coverage --cover-package=abipy --logging-level=INFO

# This is to run the integration tests (slow)
#; pytest -v abipy/integration_tests

# Generate documentation
if [[ "${PYTHON_VERSION}" == "2.7" ]]; then
    ./docs/install_reqs.sh;
    cd ./docs && make && cd ..;
fi

