"""Tests for electrons.gw module"""
from __future__ import print_function, division, absolute_import, unicode_literals

import os
import collections
import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.electrons.gw import *
from abipy.core.testing import AbipyTest


class TestQPList(AbipyTest):

    def setUp(self):
        self.sigres = sigres = abilab.abiopen(abidata.ref_file("tgw1_9o_DS4_SIGRES.nc"))
        repr(self.sigres); str(self.sigres)
        assert self.sigres.to_string(verbose=2)
        self.qplist = sigres.get_qplist(spin=0, kpoint=sigres.gwkpoints[0])

    def tearDown(self):
        self.sigres.close()

    def test_qplist(self):
        """Test QPList object."""
        qplist = self.qplist
        assert isinstance(qplist, collections.Iterable)
        self.serialize_with_pickle(qplist, protocols=[-1])

        repr(qplist); str(qplist)
        qplist_copy = qplist.copy()
        assert qplist_copy == qplist

        qpl_e0sort = qplist.sort_by_e0()
        assert qpl_e0sort.is_e0sorted
        e0mesh = qpl_e0sort.get_e0mesh()
        assert e0mesh[-1] > e0mesh[0]
        values = qpl_e0sort.get_field("qpeme0")
        assert len(values) == len(qpl_e0sort)

        qp = qpl_e0sort[2]
        value = qpl_e0sort.get_skb_field(qp.skb, "qpeme0")
        assert qp.qpeme0 == value

        with self.assertRaises(ValueError):
            qplist.get_e0mesh()
        with self.assertRaises(ValueError):
            qplist.merge(qpl_e0sort)
        with self.assertRaises(ValueError):
            qplist.merge(qplist)

        other_qplist = self.sigres.get_qplist(spin=0, kpoint=self.sigres.gwkpoints[1])
        qpl_merge = qplist.merge(other_qplist)

        assert all(qp in qpl_merge for qp in qplist)
        assert all(qp in qpl_merge for qp in other_qplist)

        # Test QPState object.
        qp = qplist[0]
        repr(qp); str(qp)
        #qp.to_string(verbose=2, title="QP State")
        assert str(qp.tips)
        assert qp.spin == 0
        assert qp.kpoint == self.sigres.gwkpoints[0]
        assert qp.kpoint is self.sigres.gwkpoints[0]

        self.assert_equal(qp.re_qpe + 1j * qp.imag_qpe, qp.qpe)
        self.assert_almost_equal(qp.e0, -5.04619941555265, decimal=5)
        self.assert_almost_equal(qp.qpe.real, -4.76022137474714)
        self.assert_almost_equal(qp.qpe.imag, -0.011501666037697)
        self.assert_almost_equal(qp.sigxme, -16.549383605401)


class TestSigresFile(AbipyTest):

    def test_readall(self):
        for path in abidata.SIGRES_NCFILES:
            with abilab.abiopen(path) as sigres:
                repr(sigres); str(sigres)
                assert sigres.to_string(verbose=2)
                assert len(sigres.structure)

    def test_base(self):
        """Test SIGRES File."""
        sigres = abilab.abiopen(abidata.ref_file("tgw1_9o_DS4_SIGRES.nc"))
        assert sigres.nsppol == 1
        sigres.print_qps(precision=5, ignore_imag=False)
        assert sigres.params["nsppol"] == sigres.nsppol
        assert not sigres.has_spectral_function

        # In this run IBZ = kptgw
        assert len(sigres.ibz) == 6
        assert sigres.gwkpoints == sigres.ibz
        # No spectral function
        assert not sigres.reader.has_spfunc
        with self.assertRaises(ValueError):
            sigres.read_sigee_skb(0, 0, 0)

        kptgw_coords = np.reshape([
            -0.25, -0.25, 0,
            -0.25, 0.25, 0,
            0.5, 0.5, 0,
            -0.25, 0.5, 0.25,
            0.5, 0, 0,
            0, 0, 0
        ], (-1, 3))

        self.assert_almost_equal(sigres.ibz.frac_coords, kptgw_coords)

        qpgaps = [3.53719151871085, 4.35685250045637, 4.11717896881632,
                  8.71122659251508, 3.29693118466282, 3.125545059031]
        self.assert_almost_equal(sigres.qpgaps, np.reshape(qpgaps, (1, 6)))

        ik = 2
        df = sigres.get_dataframe_sk(spin=0, kpoint=ik)
        same_df = sigres.get_dataframe_sk(spin=0, kpoint=sigres.gwkpoints[ik])

        assert np.all(df["qpe"] == same_df["qpe"])

        # Ignore imaginary part.
        df_real = sigres.get_dataframe_sk(spin=0, kpoint=ik, ignore_imag=True)
        assert np.all(df["qpe"].real == df_real["qpe"])

        full_df = sigres.to_dataframe()

        marker = sigres.get_marker("qpeme0")
        assert marker and len(marker.x)

        if self.has_matplotlib():
            assert sigres.plot_qps_vs_e0(fontsize=8, e0=1.0, xlims=(-10, 10), show=False)
            with self.assertRaises(ValueError):
                sigres.plot_qps_vs_e0(with_fields="qqeme0", show=False)
            assert sigres.plot_qps_vs_e0(with_fields="qpeme0", show=False)
            assert sigres.plot_qps_vs_e0(exclude_fields=["vUme"], show=False)
            assert sigres.plot_ksbands_with_qpmarkers(qpattr="sigxme", e0=None, fact=1000, show=False)
            assert sigres.plot_qpbands_ibz(show=False)

            assert sigres.plot_eigvec_qp(spin=0, kpoint=0, show=False)
            assert sigres.plot_eigvec_qp(spin=0, kpoint=None, show=False)

            assert sigres.plot_qpgaps(plot_qpmks=True, show=False)
            assert sigres.plot_qpgaps(plot_qpmks=False, show=False)

        if self.has_nbformat():
            sigres.write_notebook(nbpath=self.get_tmpname(text=True))

        sigres.close()

    def test_sigres_with_spectral_function(self):
        """Test methods to plot spectral function from SIGRES."""
        filepath = abidata.ref_file("al_g0w0_sigmaw_SIGRES.nc")
        with abilab.abiopen(filepath) as sigres:
            assert sigres.reader.has_spfunc
            sigma = sigres.read_sigee_skb(spin=0, kpoint=(0, 0, 0), band=0)
            repr(sigma); str(sigma)
            assert sigma.to_string(verbose=2)

            if self.has_matplotlib():
                assert sigma.plot(what_list="spfunc", xlims=(-10, 10), fontsize=12, show=False)
                assert sigres.plot_spectral_functions(show=False)
                assert sigres.plot_spectral_functions(include_bands=range(0, 4), show=False)

            with abilab.SigresRobot() as robot:
                robot.add_file("foo", filepath)
                robot.add_file("same", filepath)
                if self.has_matplotlib():
                    assert robot.plot_selfenergy_conv(0, (0.5, 0, 0), band=3, show=False)
                    assert robot.plot_selfenergy_conv(0, (0.5, 0, 0), band=3, sortby="nkpt", hue="nband", show=False)
                    with self.assertRaises(AttributeError):
                        assert robot.plot_selfenergy_conv(0, (0.5, 0, 0), band=3, sortby="foonkpt", hue="nband", show=False)

    def test_interpolator(self):
        """Test QP interpolation."""
        # Get quasiparticle results from the SIGRES.nc database.
        sigres = abilab.abiopen(abidata.ref_file("si_g0w0ppm_nband30_SIGRES.nc"))

        # Interpolate QP corrections and apply them on top of the KS band structures.
        # QP band energies are returned in r.qp_ebands_kpath and r.qp_ebands_kmesh.

        # Just to test call without ks_ebands.
        r = sigres.interpolate(lpratio=5,
                               ks_ebands_kpath=None,
                               ks_ebands_kmesh=None,
                               verbose=0, filter_params=[1.0, 1.0], line_density=10)

        r = sigres.interpolate(lpratio=5,
                               ks_ebands_kpath=abidata.ref_file("si_nscf_GSR.nc"),
                               ks_ebands_kmesh=abidata.ref_file("si_scf_GSR.nc"),
                               verbose=0, filter_params=[1.0, 1.0], line_density=10)

        assert r.qp_ebands_kpath is not None
        assert r.qp_ebands_kpath.kpoints.is_path
        #print(r.qp_ebands_kpath.kpoints.ksampling, r.qp_ebands_kpath.kpoints.mpdivs_shifts)
        assert r.qp_ebands_kpath.kpoints.mpdivs_shifts == (None, None)

        assert r.qp_ebands_kmesh is not None
        assert r.qp_ebands_kmesh.kpoints.is_ibz
        assert r.qp_ebands_kmesh.kpoints.ksampling is not None
        assert r.qp_ebands_kmesh.kpoints.is_mpmesh
        qp_mpdivs, qp_shifts = r.qp_ebands_kmesh.kpoints.mpdivs_shifts
        assert qp_mpdivs is not None
        assert qp_shifts is not None
        ks_mpdivs, ks_shifts = r.ks_ebands_kmesh.kpoints.mpdivs_shifts
        self.assert_equal(qp_mpdivs, ks_mpdivs)
        self.assert_equal(qp_shifts, ks_shifts)

        # Get DOS from interpolated energies.
        ks_edos = r.ks_ebands_kmesh.get_edos()
        qp_edos = r.qp_ebands_kmesh.get_edos()

        r.qp_ebands_kmesh.to_bxsf(self.get_tmpname(text=True))

        points = sigres.get_points_from_ebands(r.qp_ebands_kpath, size=24, verbose=2)

        # Plot the LDA and the QPState band structure with matplotlib.
        plotter = abilab.ElectronBandsPlotter()
        plotter.add_ebands("LDA", r.ks_ebands_kpath, edos=ks_edos)
        plotter.add_ebands("GW (interpolated)", r.qp_ebands_kpath, edos=qp_edos)

        if self.has_matplotlib():
            assert plotter.combiplot(title="Silicon band structure", show=False)
            assert plotter.gridplot(title="Silicon band structure", show=False)
            assert r.qp_ebands_kpath.plot(points=points, show=False)

        sigres.close()


class SigresRobotTest(AbipyTest):

    def test_sigres_robot(self):
        """Testing SIGRES robot."""
        filepaths = abidata.ref_files(
            "si_g0w0ppm_nband10_SIGRES.nc",
            "si_g0w0ppm_nband20_SIGRES.nc",
            "si_g0w0ppm_nband30_SIGRES.nc",
        )
        assert abilab.SigresRobot.class_handles_filename(filepaths[0])
        assert len(filepaths) == 3

        with abilab.SigresRobot.from_files(filepaths, abspath=True) as robot:
            assert robot.start is None
            start = robot.trim_paths(start=None)
            assert robot.start == start
            for p, _ in robot.items():
                assert p == os.path.relpath(p, start=start)

            assert robot.EXT == "SIGRES"
            repr(robot); str(robot)
            assert robot.to_string(verbose=2)
            assert robot._repr_html_()

            df_params = robot.get_params_dataframe()
            self.assert_equal(df_params["nsppol"].values, 1)

            label_ncfile_param = robot.sortby("nband")
            assert [t[2] for t in label_ncfile_param] == [10, 20, 30]
            label_ncfile_param = robot.sortby(lambda ncfile: ncfile.ebands.nband, reverse=True)
            assert [t[2] for t in label_ncfile_param] == [30, 20, 10]

            df_sk = robot.merge_dataframes_sk(spin=0, kpoint=[0, 0, 0])
            qpdata = robot.get_qpgaps_dataframe(with_geo=True)

            # Test plotting methods.
            if self.has_matplotlib():
                assert robot.plot_qpgaps_convergence(plot_qpmks=False, sortby=None, hue=None, show=False)
                assert robot.plot_qpgaps_convergence(plot_qpmks=True, sortby="nband", hue="ecuteps", show=False)

                assert robot.plot_qpdata_conv_skb(spin=0, kpoint=(0, 0, 0), band=3, show=False)
                assert robot.plot_qpdata_conv_skb(spin=0, kpoint=(0, 0, 0), band=5,
                        sortby="sigma_nband", hue="ecuteps", show=False)
                with self.assertRaises(TypeError):
                    robot.plot_qpdata_conv_skb(spin=0, kpoint=(0, 0, 0), band=5,
                            sortby="sigma_nband", hue="fooecueps", show=False)

                # Test plot_qpfield_vs_e0
                assert robot.plot_qpfield_vs_e0("qpeme0", sortby=None, hue=None, e0="fermie",
                        colormap="viridis", show=False)
                assert robot.plot_qpfield_vs_e0("ze0", itemp=1, sortby="ebands.nkpt", hue="scr_nband",
                        colormap="viridis", show=False)

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))

            robot.pop_label(os.path.relpath(filepaths[0], start=start))
            assert len(robot) == 2
            robot.pop_label("foobar")
            new2old = robot.change_labels(["hello", "world"], dryrun=True)
            assert len(new2old) == 2 and "hello" in new2old

            new2old = robot.remap_labels(lambda af: af.filepath, dryrun=False)
            assert len(new2old) == 2
            assert all(key == abifile.filepath for key, abifile in robot.items())
