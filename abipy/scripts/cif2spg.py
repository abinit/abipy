#!/usr/bin/env python
"""
FIXME Add doc
"""
from __future__ import print_function, division, unicode_literals

import os
import argparse
import abipy.abilab as abilab
from pymatgen.symmetry.finder import SymmetryFinder


def main():
    def str_examples():
        examples = """
    Usage example:\n
        cif2spg.py silicon.cif     => Open CIF file and visualize info on the spacegroup.
        cif2spg.py -g silicon.cif  => Same as above but use the GUI. 
    """
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: 
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")

        sys.exit(error_code)


    parser = argparse.ArgumentParser(epilog=str_examples(), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('cif_file', nargs=1, help="CIF File")

    parser.add_argument('-g', '--gui', action="store_true", help="Enable GUI")

    # Parse command line.
    try:
        options = parser.parse_args()
    except: 
        show_examples_and_exit(error_code=1)

    cif_file = options.cif_file[0]
    structure = abilab.Structure.from_file(cif_file)
    #print(structure.to_abivars())

    finder = SymmetryFinder(structure, symprec=1e-5, angle_tolerance=5)

    data = finder.get_symmetry_dataset()
    from pprint import pprint

    if not options.gui:
        pprint(data)

    else:
        from StringIO import StringIO
        import wx
        from abipy.gui.editor import SimpleTextViewer

        stream = StringIO()
        pprint(data, stream=stream)
        stream.seek(0)
        text = "".join(stream)
        app = wx.App()
        frame = SimpleTextViewer(None, text=text)
        frame.Show()
        app.MainLoop()

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
