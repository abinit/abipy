#!/bin/bash
set -ev  # exit on first error, print each command

echo "PMG_MAPI_KEY: 8pkvwRLQSCVbW2Fe" > ${HOME}/.pmgrc.yaml

abicheck.py --with-flow

nosetests -v --with-coverage --cover-package=abipy --logging-level=INFO
#nosetests abipy -v --with-coverage --cover-package=abipy --logging-level=INFO

#pytest --cov-config .coveragerc --cov=abipy -v  abipy # --doctest-modules 
# This is to run the integration tests (append results)
# pytest --cov-config .coveragerc --cov=abipy --cov-append -v abipy/integration_tests

# Generate documentation
if [[ "${PYTHON_VERSION}" == "2.7" && "${TRAVIS_OS_NAME}" == "linux" ]]; then
    ./docs/install_reqs.sh;
    cd ./docs && make && cd ..;
fi
