"""Tests for electrons.bse module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.electrons.scissors import *
from abipy import abilab


class TestScissors(AbipyTest):

    def test_scissors_polyfit(self):
        """Testing scissors from SIGRES file."""
        # Get the quasiparticle results from the SIGRES.nc database.
        with abilab.abiopen(abidata.ref_file("si_g0w0ppm_nband30_SIGRES.nc")) as sigma_file:
            qplist_spin = sigma_file.qplist_spin

        # Construct the scissors operator
        domains = [[-10, 6.1], [6.1, 18]]
        scissors = qplist_spin[0].build_scissors(domains, bounds=None)
        #scissors = qplist_spin[0].build_scissors(domains, bounds=None, plot=True)

        # Read the KS band energies computed on the k-path
        with abilab.abiopen(abidata.ref_file("si_nscf_GSR.nc")) as nc:
           ks_bands = nc.ebands

        # Read the KS band energies computed on the Monkhorst-Pack (MP) mesh
        # and compute the DOS with the Gaussian method
        with abilab.abiopen(abidata.ref_file("si_scf_GSR.nc")) as nc:
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
            plotter = abilab.ElectronBandsPlotter()
            plotter.add_ebands("LDA", ks_bands, edos=ks_dos)
            plotter.add_ebands("LDA+scissors(e)", qp_bands, edos=qp_dos)

            # By default, the two band energies are shifted wrt to *their* fermi level.
            # Use e=0 if you don't want to shift the eigenvalus
            # so that it's possible to visualize the QP corrections.
            plotter.combiplot(title="Silicon band structure", show=False)
            plotter.gridplot(title="Silicon band structure", show=False)

    def test_builder(self):
        """Testing ScissorsBuilder."""
        builder = ScissorsBuilder.from_file(abidata.ref_file("si_g0w0ppm_nband30_SIGRES.nc"))
        assert builder.nsppol == 1

        # To select the domains esplicitly (optional)
        builder.build(domains_spin=[[-10, 6.02], [6.1, 20]])

        # Test pickle.
        tmp_path = self.get_tmpname(suffix=".pickle")
        builder.pickle_dump(tmp_path)
        new = ScissorsBuilder.pickle_load(tmp_path)

        if self.has_matplotlib():
            # To plot the QP results as function of the KS energy:
            assert builder.plot_qpe_vs_e0(show=False)
            # To compare the fitted results with the ab-initio data:
            assert builder.plot_fit(show=False)

            # To plot the corrected bands:
            bands_filepath = abidata.ref_file("si_nscf_WFK.nc")
            assert builder.plot_qpbands(bands_filepath, show=False)

            dos_filepath = abidata.ref_file("si_scf_GSR.nc")
            assert builder.plot_qpbands(bands_filepath, dos_filepath=dos_filepath,
                    dos_args=None, show=False)
