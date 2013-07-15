#!/usr/bin/env python
"""
This scripts runs all the python scripts located in this directory allowing
the user to change the matplotlib backend. 

Usage: 
    run.py [backend]
"""
from __future__ import print_function

import sys
import os 
import shutil

CONF_FILE = None
BKP_FILE = None

def change_backend(new_backend=None):
    """Change the backend by modifying the matplotlib configuration file."""
    global CONF_FILE, BKP_FILE

    if new_backend is None:
        return 

    home = os.environ["HOME"]
    CONF_FILE = conf_file = os.path.join(home, ".matplotlib", "matplotlibrc")

    BKP_FILE = conf_file + ".bkp"

    if os.path.exists(conf_file):
        shutil.copy(conf_file, BKP_FILE)

        with open(conf_file, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if line.strip().startswith("backend"):
                lines[i] = "backend : " + new_backend + "\n"

        with open(conf_file, "w") as f:
            f.writelines(lines)


def revert_backend():
    global CONF_FILE, BKP_FILE
    #print("reverting: BKP_FILE %s --> CONF %s" % (BKP_FILE, CONF_FILE))
    if BKP_FILE is not None:
        shutil.move(BKP_FILE, CONF_FILE)


def main():
    new_backend = None
    if len(sys.argv) > 1:
        new_backend = sys.argv[1]

    change_backend(new_backend=new_backend)

    import subprocess
    for fname in os.listdir("."):
        if fname.endswith(".py") and fname != "run.py":
            retcode = subprocess.call(["python", fname])
            if retcode != 0: break

    revert_backend()
    return retcode

if __name__ == "__main__":
    sys.exit(main())
