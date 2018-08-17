#!/usr/bin/env python
r"""
Phonopy + Abinit Flow
=====================

This example shows how to compute phonon frequencies with phonopy (supercells and finite-difference method).
This approach could be useful to obtain vibrational properties with XC functionals for which DFPT is not yet implemented.

.. warning:

    This example requires the `phonopy package <http://atztogo.github.io/phonopy/examples.html>`_
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.flowtk.abiphonopy import PhonopyWork


def build_flow(options):
    """
    Create a `Flow` for phonon calculations with phonopy:
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize structure and pseudos
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
    pseudos = abidata.pseudos("14si.pspnc")

    # Build input for GS calculation.
    gsinp = abilab.AbinitInput(structure, pseudos)
    gsinp.set_vars(ecut=4, nband=4, toldff=1.e-6)

    # This gives ngkpt = 4x4x4 with 4 shifts for the initial unit cell.
    # The k-point sampling will be rescaled when we build the supercell in PhonopyWork.
    gsinp.set_autokmesh(nksmall=4)
    #gsinp.set_vars(ngkpt=[4, 4, 4])

    flow = flowtk.Flow(workdir=options.workdir)

    # Use a 2x2x2 supercell to compute phonons with phonopy
    work = PhonopyWork.from_gs_input(gsinp, scdims=[2, 2, 2])
    flow.register_work(work)

    return flow


# This block generates the thumbnails in the Abipy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("READTHEDOCS", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    #build_flow(options).plot_networkx(with_edge_labels=True, tight_layout=True)
    build_flow(options).graphviz_imshow()


@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())

############################################################################
#
# Run the script with:
#
#     run_phonopy_si.py -s
#
# the output results are produced in ``flow_phonopy_si/w0/outdata/``
# Follow the instructions in the README file:
#
# .. code-block:: md
#
#    To plot bands, use:
#            phonopy -p band.conf
#
#    To plot phonon dos, use:
#            phonopy -p dos.conf
#
#    To plot bands and dos, use:
#            phonopy -p band-dos.conf
#
#    See also:
#            http://atztogo.github.io/phonopy/examples.html
#            http://atztogo.github.io/phonopy/setting-tags.html#setting-tags
#
# The command:
#
#   phonopy -p band-dos.conf
#
# will produce:
#
# .. image:: https://github.com/abinit/abipy_assets/blob/master/run_phonopy_si.png?raw=true
#    :alt: Phonon Band structure computed with phonopy.
#
