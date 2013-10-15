#!/usr/bin/env python
from __future__ import division, print_function

import collections
import yaml
from pymatgen.io.abinitio.events import YamlFileReader

if __name__ == "__main__":
    import sys
    log = sys.argv[1]

    with YamlFileReader(log) as r:
        all_docs = r.all_yaml_docs()

        for doc in all_docs:
            #print(80*"*")
            #print(doc)
            #print(80*"*")
            obj = yaml.load(doc)
            print(obj)
