#!/usr/bin/env python
"""
This script generates RST files to be included in the AbiPy website.
It is automatically executed by make
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os

ABIPY_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "abipy")

def generate_manager_examples():
    rst_path = os.path.join("workflows", "manager_examples.rst")
    print("Generating RST file %s with all manager files..." % os.path.relpath(rst_path))
    dirpath = os.path.join(ABIPY_ROOT, "data", "managers")
    assert os.path.isdir(dirpath)
    manager_files = sorted([os.path.join(dirpath, f) for f in os.listdir(dirpath) if f.endswith("_manager.yml")])
    assert manager_files

    lines = []; app = lines.append
    for f in manager_files:
        machine_name = os.path.basename(f).replace("_manager.yml", "").capitalize()
        app(machine_name)
        app("-" * len(machine_name))
        app("""
.. code-block:: yaml

""")
        with open(f, "rt") as fh:
            lines.extend(["\t%s" % l.rstrip() for l in fh])
        app("\n")

    text = """\

.. _manager-examples:

****************
Manager Examples
****************

.. contents::
   :backlinks: top

""" + "\n".join(lines)

    rst_path = os.path.join("workflows", "manager_examples.rst")
    with open(rst_path, "wt") as fh:
        fh.write(text)


def main():
    generate_manager_examples()
    return 0


if __name__ == "__main__":
    sys.exit(main())
