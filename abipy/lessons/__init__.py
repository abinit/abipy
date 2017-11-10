"""Directory with abipy + abinit tutorials."""
from __future__ import division, print_function, unicode_literals, absolute_import

import os

from importlib import import_module


def help():
    """List the available tutorials with a brief description."""
    docs = [m.__doc__ for mod in get_mods()]
    return "\n".join(docs)


def get_mods():
    mods = []
    dirpath = os.path.abspath(os.path.dirname(__file__))
    for pyfile in os.listdir(dirpath):
        if not (pyfile.startswith("lesson") and pyfile.endswith(".py")): continue
        path = os.path.join(os.path.dirname(dirpath), pyfile)
        mods.append(import_module(pyfile.replace(".py", "")))

    return mods


def generate_manfiles():
    """Generate the man files from the lesson scripts."""
    for mod in get_mods():
        try:
            lesson = mod.Lesson()

        except AttributeError:
            print("Exception in mod %s" % mod.__name__)
            continue

        lesson._gen_manfile()
        print("[%s] man file generated" % mod.__name__)


#if __name__ == "__main__":
#    generate_manfiles()
