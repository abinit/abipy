#!/bin/bash
set -ev  # exit on first error, print each command

abicheck.py --with-flow

nosetests -v --with-coverage --cover-package=abipy --logging-level=INFO
#nosetests abipy -v --with-coverage --cover-package=abipy --logging-level=INFO

# This is to run the integration tests (slow)
# pytest -v --cov=abipy--doctest-modules --durations=10 --ignore=./docs/ abipy/integration_tests
# pytest abipy/integration_tests --ignore=./docs/

# Generate documentation
if [[ "${PYTHON_VERSION}" == "2.7" && "${TRAVIS_OS_NAME}" == "linux" ]]; then
    ./docs/install_reqs.sh;
    cd ./docs && make && cd ..;
fi
