#!/usr/bin/env python
"""This script fetches the official yaml files from the ABINIT source tree."""
from __future__ import print_function, division, unicode_literals

import sys
import os
import shutil

if "-h" in sys.argv or "--help" in sys.argv:
    print("usage: update ~abinit/doc/input_variables/")
    sys.exit(1)


if len(sys.argv) > 1:
    doc_dir = sys.argv[1]
else:
    doc_dir = os.path.expanduser("~/git/abinit/doc/input_variables/")
    if not os.path.exists(doc_dir):
        print("usage: update ~abinit/doc/input_variables/")
        sys.exit(1)


print("Will get abinit variables from:", doc_dir)
yaml_basenames = ["abinit_vars.yml", "characteristics.yml", "sections.yml"]

yaml_paths = [os.path.join(doc_dir, f) for f in yaml_basenames] 

for p in yaml_paths:
    if not os.path.exists(p):
        raise ValueError("Wrong ABINIT doc dir %s\n" % doc_dir + 
                         "Cannot find %s:\n" % p)

for p in yaml_paths:
    shutil.copy(p, os.path.basename(p))

# Update the json file with the list of varnames.
# Read new yaml file and dump key list in json format.
import json
import yaml
import abipy.abio.abivars_db

print("Updating name list in abinit_vars.yml ...")

with open("abinit_vars.yml", "r") as fh:
    varlist = yaml.load(fh)
    varnames = sorted([v.varname for v in varlist])

with open("abinit_vars.json", "wt") as fh:
    json.dump(varnames, fh)
