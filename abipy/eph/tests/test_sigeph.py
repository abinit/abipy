"""Tests for sigeph module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import collections
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy import abilab


class SigEPhFileTest(AbipyTest):

    def test_sigeph_file(self):
        """Tests for SigEPhFile."""
        sigeph = abilab.abiopen(abidata.ref_file("diamond_444q_SIGEPH.nc"))
        repr(sigeph); str(sigeph)
        assert sigeph.to_string(verbose=2)
        assert sigeph.nsppol == 1 and sigeph.nspden == 1 and sigeph.nspinor == 1
        assert sigeph.structure.formula == "C2"
        assert sigeph.ebands.kpoints.is_ibz
        self.assert_equal(sigeph.ebands.kpoints.ksampling.mpdivs, [4, 4, 4])

        assert sigeph.nkcalc == 2
        self.assert_equal(sigeph.ngqpt.flatten(), [4, 4, 4])
        assert not sigeph.imag_only
        assert sigeph.symsigma == 0
        assert sigeph.ntemp ==  6
        assert sigeph.nband == 54
        assert sigeph.nqbz == sigeph.ngqpt.prod()
        assert sigeph.nqibz == 8
        # FIXME
        #self.assert_almost_equal(sigeph.zcut, 0.001)
        assert sigeph.has_spectral_function and sigeph.reader.nwr == 101
        assert not sigeph.has_eliashberg_function
        assert len(sigeph.mu_e) == sigeph.ntemp
        assert "nbsum" in sigeph.params
        assert "eph_fsewin" in sigeph.params

        assert sigeph.ks_dirgaps.shape == (sigeph.nsppol, sigeph.nkcalc)
        assert sigeph.qp_dirgaps_t.shape == (sigeph.nsppol, sigeph.nkcalc, sigeph.ntemp)

        assert len(sigeph.sigma_kpoints) == 2
        assert not sigeph.sigma_kpoints.is_ibz and not sigeph.sigma_kpoints.is_path
        assert sigeph.sigma_kpoints[0] == [0, 0, 0]
        assert sigeph.sigma_kpoints[1] == [0.5, 0, 0]
        assert sigeph.bstart_sk.shape == (sigeph.nsppol, sigeph.nkcalc)
        self.assert_equal(sigeph.bstart_sk.flatten(), [0, 0])
        self.assert_equal(sigeph.nbcalc_sk.flatten(), [8, 8])
        assert sigeph.bstart_sk.shape == sigeph.nbcalc_sk.shape
        self.assert_equal(sigeph.bstop_sk.flatten(), [8, 8])
        self.assert_almost_equal(sigeph.tmesh[0], 5)
        assert sigeph.nbsum == 54
        assert sigeph.sigkpt2index(0) == 0
        assert sigeph.sigkpt2index([0.5, 0.0, 0.0]) == 1

        # Test find_qpkinds
        kpt, ikc = sigeph.find_qpkinds(1)[0]
        assert ikc == 1 and kpt == sigeph.sigma_kpoints[ikc]
        kpt, ikc = sigeph.find_qpkinds([[0.5, 0, 0]])[0]
        assert ikc == 1 and kpt == sigeph.sigma_kpoints[ikc]
        kpt_list, ikc_list = zip(*sigeph.find_qpkinds([0, 1]))
        assert len(kpt_list) == 2
        assert ikc_list == (0, 1)
        qpkinds = sigeph.find_qpkinds([[0, 0, 0], [0.5, 0, 0]])
        assert len(qpkinds) == 2
        assert tuple(t[1] for t in qpkinds) == (0, 1)
        assert sigeph.find_qpkinds(qpkinds) is qpkinds

        # Test ksampling
        ksamp = sigeph.reader.read_ksampling_info()
        assert ksamp.is_mesh and not ksamp.is_path
        assert ksamp.has_diagonal_kptrlatt
        self.assert_equal(ksamp.mpdivs, [4, 4, 4])
        self.assert_equal(ksamp.shifts.ravel(), [0, 0, 0])
        assert ksamp.to_string(title="Ksampling")

        # Test Dataframe construction.
        data_sk = sigeph.get_dataframe_sk(spin=0, kpoint=[0.5, 0.0, 0.0], with_spin=True)
        assert "qpeme0" in data_sk
        assert np.all(data_sk["spin"] == 0)
        #self.assert_almost_equal(data_sk["kpoint"].values[0].frac_coords, [0.5, 0.0, 0.0])

        itemp = 0
        data_sk = sigeph.get_dataframe_sk(spin=0, kpoint=[0.5, 0.0, 0.0], itemp=itemp)
        self.assert_equal(data_sk["tmesh"].values, sigeph.tmesh[itemp])

        data = sigeph.get_dataframe()
        assert "ze0" in data

        if self.has_matplotlib():
            # Test sigeph plot methods.
            assert sigeph.plot_qpgaps_t(show=False)
            assert sigeph.plot_qpgaps_t(plot_qpmks=True, show=False)
            assert sigeph.plot_qpdata_t(spin=0, kpoint=(0, 0, 0), show=False)
            assert sigeph.plot_qpbands_ibzt(itemp_list=[0, 1], colormap="viridis", show=False)

        if self.has_nbformat():
            sigeph.write_notebook(nbpath=self.get_tmpname(text=True))

        # Test self-energy object.
        with self.assertRaises(ValueError):
            sigeph.reader.sigkpt2index([0.3, 0.5, 0.4])
        with self.assertRaises(ValueError):
            sigeph.reader.read_sigeph_skb(spin=0, kpoint=5, band=0)
        with self.assertRaises(ValueError):
            sigeph.reader.read_sigeph_skb(spin=0, kpoint=[0.3, 0.5, 0.4], band=0)
        with self.assertRaises(ValueError):
            sigeph.reader.read_sigeph_skb(spin=0, kpoint=0, band=100)

        sigma = sigeph.reader.read_sigeph_skb(spin=0, kpoint=[0.5, 0, 0], band=3)
        repr(sigma); str(sigma)
        assert sigma.to_string(verbose=2, title="Fan-Migdal Self-energy.")
        assert sigma.spin == 0 and sigma.band == 3 and sigma.kpoint == [0.5, 0, 0]
        assert sigma.tmesh is sigeph.tmesh
        assert sigma.wmesh.shape == (sigma.nwr,)
        if self.has_matplotlib():
            assert sigma.plot_tdep(show=False)

        # Test QpTempState
        qp = sigeph.reader.read_qp(spin=0, kpoint=0, band=3, ignore_imag=False)
        repr(qp); str(qp)
        assert qp.to_string(verbose=2, title="QPTempState")
        assert qp._repr_html_()
        assert qp.spin == 0 and qp.kpoint == [0, 0, 0] and qp.band == 3
        assert qp.skb == (0, [0, 0, 0], 3)
        #assert len(qp.qpeme0) == qpt.ntemp
        self.assert_equal(qp.qpeme0, (qp.qpe - qp.e0).real)
        self.assert_equal(qp.re_qpe + 1j * qp.imag_qpe, qp.qpe)
        self.assert_equal(qp.re_fan0 + 1j * qp.imag_fan0, qp.fan0)
        fields = qp.get_fields()
        assert all(k in fields for k in ("qpeme0", "kpoint"))
        #self.assert_almost_equal(qp.e0, -5.04619941555265, decimal=5)
        #self.assert_almost_equal(qp.qpe.real, -4.76022137474714)
        #self.assert_almost_equal(qp.qpe.imag, -0.011501666037697)
        #self.assert_almost_equal(qp.sigxme, -16.549383605401)
        if self.has_matplotlib():
            # Test plot methods.
            assert qp.plot(show=False)

        # Test QPList
        qplist = sigeph.reader.read_qplist_sk(spin=0, kpoint=[0, 0, 0], ignore_imag=False)
        assert isinstance(qplist, collections.Iterable)
        # TODO
        #self.serialize_with_pickle(qplist, protocols=[-1])

        repr(qplist); str(qplist)
        assert qplist.to_string(verbose=2, title="QP list")

        #qplist_copy = qplist.copy()
        #assert qplist_copy == qplist
        qpl_e0sort = qplist.sort_by_e0()
        assert qpl_e0sort.is_e0sorted

        e0mesh = qpl_e0sort.get_e0mesh()
        assert e0mesh[-1] > e0mesh[0]
        #values = qpl_e0sort.get_field("qpeme0")
        #assert len(values) == len(qpl_e0sort)

        #qp = qpl_e0sort[2]
        #value = qpl_e0sort.get_skb_field(qp.skb, "qpeme0")
        #assert qp.qpeme0 == value

        with self.assertRaises(ValueError):
            qplist.get_e0mesh()
        with self.assertRaises(ValueError):
            qplist.merge(qpl_e0sort)
        with self.assertRaises(ValueError):
            qplist.merge(qplist)

        other_qplist = sigeph.reader.read_qplist_sk(spin=0, kpoint=sigeph.sigma_kpoints[1])
        qpl_merge = qplist.merge(other_qplist)

        assert all(qp in qpl_merge for qp in qplist)
        assert all(qp in qpl_merge for qp in other_qplist)

        if self.has_matplotlib():
            assert qpl_merge.plot_vs_e0(show=False)

        qplist_spin = sigeph.qplist_spin
        assert len(qplist_spin) == sigeph.nsppol
        if self.has_matplotlib():
            assert sigeph.plot_qps_vs_e0(show=False)

        print(sigeph.kcalc2ibz)
        dos = sigeph.get_linewidth_dos()

        sigeph.close()

    # TODO: Need new files with IBZ.
    def test_sigeph_interpolation(self):
        """Test interpolation and TdepElectronBands"""
        sigeph = abilab.abiopen(abidata.ref_file("diamond_444q_full_SIGEPH.nc"))

        # Test interpolation without KS bands.
        tdep_nopath = sigeph.interpolate(itemp_list=0)
        repr(tdep_nopath); str(tdep_nopath)
        assert tdep_nopath.to_string(verbose=2)
        assert tdep_nopath.ntemp == 1
        #assert not tdep_nopath.has_kpath
        #assert not tdep_nopath.ks_ebands_kpath is None
        #assert not tdep_nopath.has_kmesh
        #assert not tdep_nopath.ks_ebands_kmesh is None
        same_tdep_nopath = tdep_nopath.__class__.pickle_load(tdep_nopath.pickle_dump())

        if self.has_matplotlib():
            assert tdep_nopath.plot_itemp(itemp=0, fontsize=8, show=False)
            assert tdep_nopath.plot_itemp_with_lws_vs_e0(itemp=0, fontsize=8, with_ratios=(3, 1), show=False)
            assert tdep_nopath.plot(show=False)
            assert tdep_nopath.plot_lws_vs_e0(itemp_list=[0, -1], show=False)
            assert tdep_nopath.get_ebands_plotter()
            #assert tdep_nopath.get_edos_plotter() is None

        # Test interpolation with KS bands.
        tdep = sigeph.interpolate(itemp_list=None, lpratio=5, ks_ebands_kpath=None, ks_ebands_kmesh=None, ks_degatol=1e-4,
                    vertices_names=None, line_density=20, filter_params=None, only_corrections=False, verbose=0)
        repr(tdep); str(tdep)
        assert tdep.to_string(verbose=2)
        assert tdep_nopath.ntemp == 1
        #assert not tdep_nopath.ks_ebands_kpath is None
        #assert tdep.has_kpath
        #assert not tdep_nopath.has_kmesh
        #assert not tdep_nopath.ks_ebands_kmesh is None

        if self.has_matplotlib():
            assert tdep_nopath.plot_itemp(itemp=0, fontsize=8, show=False)
            assert tdep_nopath.plot_itemp_with_lws_vs_e0(itemp=0, fontsize=8, show=False)
            assert tdep_nopath.plot(show=False)
            assert tdep_nopath.plot_lws_vs_e0(itemp_list=[0, -1], show=False)
            assert tdep.get_ebands_plotter()
            #assert tdep.get_edos_plotter() is None

        sigeph.close()

    # TODO: Need new files with IBZ.
    def test_sigeph_boltztrap(self):
        """Test boltztrap interpolation"""
        sigeph = abilab.abiopen(abidata.ref_file("diamond_444q_full_SIGEPH.nc"))
        sigeph.get_lifetimes_boltztrap("diamond", workdir=self.mkdtemp())
        sigeph.close()

    def test_sigeph_robot(self):
        """Tests for SigEPhRobot."""
        filepaths = [
            abidata.ref_file("diamond_444q_SIGEPH.nc"),
        ]
        with abilab.SigEPhRobot.from_files(filepaths) as robot:
            robot.add_file("same_file", filepaths[0])
            repr(robot); str(robot)
            assert robot.to_string(verbose=2)
            assert len(robot) == 2

            data = robot.get_dataframe()
            assert "re_qpe" in data

            # Test plot methods
            if self.has_matplotlib():
                assert robot.plot_selfenergy_conv(spin=0, kpoint=0, band=0, show=False)
                assert robot.plot_selfenergy_conv(spin=0, kpoint=0, band=0, sortby="nbsum", hue="nqibz", show=False)
                assert robot.plot_qpgaps_t(show=False)
                assert robot.plot_qpgaps_t(plot_qpmks=True, show=False)

                assert robot.plot_qpgaps_convergence(itemp=0, sortby="nbsum", show=False)
                assert robot.plot_qpgaps_convergence(itemp=0, sortby="nbsum", hue="nqibz", show=False)

                assert robot.plot_qpdata_conv_skb(spin=0, kpoint=(0, 0, 0), band=3,
                        itemp=-1, show=False)
                assert robot.plot_qpdata_conv_skb(spin=0, kpoint=(0, 0, 0), band=3,
                        itemp=0, sortby="nbsum", hue="nqibz", show=False)

                # Test plot_qpfield_vs_e0
                assert robot.plot_qpfield_vs_e0("qpeme0", itemp=1, sortby=None, hue=None,
                        colormap="viridis", show=False)
                assert robot.plot_qpfield_vs_e0("ze0", itemp=1, sortby="nbsum", hue="nqibz",
                        colormap="viridis", show=False)

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))
