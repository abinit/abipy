#!/bin/bash
set -ev  # exit on first error, print each command

pip install -r requirements.txt
pip install -r requirements-optional.txt
pip install -r requirements-panel.txt
pip install -r requirements-tests.txt
python setup.py develop
# This to bypass a breaking API change in the buildbot testfarm
#pip install pymatgen==v2024.2.8 -U
