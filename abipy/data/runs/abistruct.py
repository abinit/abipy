#!/usr/bin/env python
from __future__ import division, print_function

import sys
import os

from abipy import abilab

def main():
    import sys
    filepath = sys.argv[1]

    structure = abilab.Structure.from_file(filepath)

    for format in ["cif", "POSCAR", "cssr", "json"]:
        s = structure.convert(format=format)

        print((" Abinit --> %s " % format).center(80, "*"))
        print(s)


if __name__ == "__main__":
    sys.exit(main())
