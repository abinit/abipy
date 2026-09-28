#!/bin/bash
set -ev  # exit on first error, print each command

# Install AbiPy and all development extras in one resolve. Dependency groups
# are defined in pyproject.toml, the single source of package metadata.
python -m pip install --editable ".[optional,panel,tests,dev]"

invoke submodules
