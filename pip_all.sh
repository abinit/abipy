#!/bin/bash
set -ev  # exit on first error, print each command

pip install -r requirements.txt
pip install -r requirements-optional.txt
pip install -r requirements-panel.txt
pip install -r requirements-tests.txt
python setup.py develop

pip install invoke
invoke submodules
