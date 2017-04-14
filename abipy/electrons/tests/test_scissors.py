"""Tests for electrons.bse module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import abipy.data as abidata

from abipy.electrons.scissors import *
from abipy.core.testing import AbipyTest

from abipy.abilab import abiopen, ElectronBandsPlotter


class TestScissors(AbipyTest):

    def test_scissors_polyfit(self):
        """Testing scissors from SIGRES file."""
        # Get the quasiparticle results from the SIGRES.nc database.
        sigma_file = abiopen(abidata.ref_file("si_g0w0ppm_nband30_SIGRES.nc"))
        qplist_spin = sigma_file.qplist_spin
        sigma_file.close()

        # Construct the scissors operator
        domains = [[-10, 6.1], [6.1, 18]]
        scissors = qplist_spin[0].build_scissors(domains, bounds=None)
        #scissors = qplist_spin[0].build_scissors(domains, bounds=None, plot=True)

        # Read the KS band energies computed on the k-path
        with abiopen(abidata.ref_file("si_nscf_GSR.nc")) as nc:
           ks_bands = nc.ebands

        # Read the KS band energies computed on the Monkhorst-Pack (MP) mesh
        # and compute the DOS with the Gaussian method
        with abiopen(abidata.ref_file("si_scf_GSR.nc")) as nc:
            ks_mpbands = nc.ebands
        ks_dos = ks_mpbands.get_edos()

        # Apply the scissors operator first on the KS band structure
        # along the k-path then on the energies computed with the MP mesh.
        qp_bands = ks_bands.apply_scissors(scissors)
        qp_mpbands = ks_mpbands.apply_scissors(scissors)

        # Compute the DOS with the modified QPState energies.
        qp_dos = qp_mpbands.get_edos()

        # Plot the LDA and the QPState band structure with matplotlib.
        if self.has_matplotlib():
            plotter = ElectronBandsPlotter()
            plotter.add_ebands("LDA", ks_bands, dos=ks_dos)
            plotter.add_ebands("LDA+scissors(e)", qp_bands, dos=qp_dos)

            # By default, the two band energies are shifted wrt to *their* fermi level.
            # Use e=0 if you don't want to shift the eigenvalus
            # so that it's possible to visualize the QP corrections.
            plotter.combiplot(title="Silicon band structure", show=False)
            plotter.gridplot(title="Silicon band structure", show=False)
