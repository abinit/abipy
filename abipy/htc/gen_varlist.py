#!/usr/bin/env python
"""
Script to generate the json files with the list of Abinit variables.
The names are taken from the YAML file extracted from the official documentation.
"""
from abipy.htc.abivars_db import get_abinit_variables

# Get string.
s = get_abinit_variables().json_dumps_varnames()

# Write JSON file.
with open("abinit_vars.json", "w") as fh:
    fh.write(s)
