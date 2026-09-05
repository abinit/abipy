#!/bin/bash
set -ev  # exit on first error, print each command

# Install all requirements files in ONE pip resolve. Separate `pip install -r`
# calls each resolve independently, so a package pulled in transitively by an
# earlier file (e.g. jinja2 via panel, from requirements-panel.txt) can be left
# behind a later file's own dependency's floor on that same package (e.g.
# pandas' optional jinja2>=3.1.5 requirement for pandas.io.formats.style) --
# pip never gets a chance to reconcile the two in a single cross-file resolve.
pip install -r requirements.txt -r requirements-optional.txt -r requirements-panel.txt -r requirements-tests.txt
python setup.py develop

pip install invoke
invoke submodules
