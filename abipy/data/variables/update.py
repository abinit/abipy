#!/usr/bin/env python
"""This script fetches the official yaml files from the ABINIT source tree."""

import sys
import os
import shutil

try
    doc_dir = sys.argv[1]
except IndexError:
    raise 

yaml_basenames = ["abinit_vars.yml", "characteristics.yml", "sections.yml"]

yaml_paths = [os.path.join(doc_dir, f) for f in yaml_basenames] 

for p in yaml_paths:
    if not os.path.exists(p):
        raise ValueError("Wrong ABINIT doc dir %s\n" % doc_dir + 
                         "Cannot find %s:\n" % p)

for p in yaml_paths:
    shutil.copy(p, os.path.basename(p))
