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
        #raise self.SkipTest("Disabled")
        ncfile = abilab.abiopen(abidata.ref_file("diamond_444q_SIGEPH.nc"))
        repr(ncfile); str(ncfile)
        assert ncfile.to_string(verbose=2)
        assert ncfile.nsppol == 1 and ncfile.nspden == 1 and ncfile.nspinor == 1
        assert ncfile.structure.formula == "C2"
        assert ncfile.ebands.kpoints.is_ibz
        self.assert_equal(ncfile.ebands.kpoints.ksampling.mpdivs, [4, 4, 4])

        assert ncfile.nkcalc == 2
        self.assert_equal(ncfile.ngqpt.flatten(), [4, 4, 4])
        assert ncfile.symsigma == 0
        assert ncfile.ntemp ==  6
        assert ncfile.nband == 54
        assert ncfile.nqbz == ncfile.ngqpt.prod()
        assert ncfile.nqibz == 8
        #self.assert_almost_equal(ncfile.eta, 0.001)
        assert ncfile.has_spectral_function and ncfile.reader.nwr == 101
        assert len(ncfile.mu_e) == ncfile.ntemp
        assert "nbsum" in ncfile.params

        assert ncfile.ks_dirgaps.shape == (ncfile.nsppol, ncfile.nkcalc)
        assert ncfile.qp_dirgaps_t.shape == (ncfile.nsppol, ncfile.nkcalc, ncfile.ntemp)

        assert len(ncfile.sigma_kpoints) == 2
        assert not ncfile.sigma_kpoints.is_ibz and not ncfile.sigma_kpoints.is_path
        assert ncfile.sigma_kpoints[0] == [0, 0, 0]
        assert ncfile.sigma_kpoints[1] == [0.5, 0, 0]
        assert ncfile.bstart_sk.shape == (ncfile.nsppol, ncfile.nkcalc)
        self.assert_equal(ncfile.bstart_sk.flatten(), [0, 0])
        self.assert_equal(ncfile.nbcalc_sk.flatten(), [8, 8])
        assert ncfile.bstart_sk.shape == ncfile.nbcalc_sk.shape
        self.assert_equal(ncfile.bstop_sk.flatten(), [8, 8])
        self.assert_almost_equal(ncfile.tmesh[0], 5)
        # FIXME
        #self.assert_almost_equal(ncfile.zcut, )
        assert ncfile.nbsum == 54
        assert ncfile.sigkpt2index(0) == 0
        assert ncfile.sigkpt2index([0.5, 0.0, 0.0]) == 1

        # Test Dataframe construction.
        data_sk = ncfile.get_dataframe_sk(spin=0, sigma_kpoint=[0.5, 0.0, 0.0])
        assert "qpeme0" in data_sk
        assert np.all(data_sk["spin"] == 0)
        self.assert_almost_equal(data_sk["kpoint"].values[0].frac_coords, [0.5, 0.0, 0.0])

        data = ncfile.get_dataframe()
        assert "ze0" in data

        if self.has_matplotlib():
            # Test ncfile plot methods.
            assert ncfile.plot_qpgaps_t(show=False)
            assert ncfile.plot_qpgaps_t(plot_qpmks=True, show=False)
            assert ncfile.plot_qpdata_t(spin=0, sigma_kpoint=(0, 0, 0), show=False)

        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

        # Test self-energy object.
        with self.assertRaises(ValueError):
            ncfile.reader.sigkpt2index([0.3, 0.5, 0.4])
        with self.assertRaises(ValueError):
            ncfile.reader.read_sigma_eph(spin=0, sigma_kpoint=5, band=0)
        with self.assertRaises(ValueError):
            ncfile.reader.read_sigma_eph(spin=0, sigma_kpoint=[0.3, 0.5, 0.4], band=0)
        with self.assertRaises(ValueError):
            ncfile.reader.read_sigma_eph(spin=0, sigma_kpoint=0, band=100)

        sigma = ncfile.reader.read_sigma_eph(spin=0, sigma_kpoint=[0.5, 0, 0], band=3)
        repr(sigma); str(sigma)
        assert sigma.to_string(verbose=2, title="Fan-Migdal Self-energy.")
        assert sigma.spin == 0 and sigma.band == 3 and sigma.kpoint == [0.5, 0, 0]
        assert sigma.tmesh is ncfile.tmesh
        assert sigma.wmesh.shape == (sigma.nwr,)
        if self.has_matplotlib():
            assert sigma.plot_tdep(show=False)

        # Test QpTempState
        qp = ncfile.reader.read_qp(spin=0, sigma_kpoint=0, band=3, ignore_imag=False)
        repr(qp); str(qp)
        assert qp.to_string(verbose=2, title="QPTempState")
        assert qp._repr_html_()
        assert qp.spin == 0 and qp.kpoint == [0, 0, 0] and qp.band == 3
        assert qp.skb == (0, [0, 0, 0], 3)
        #assert len(qp.qpeme0) == qpt.ntemp
        self.assert_equal(qp.qpeme0, qp.qpe - qp.e0)
        self.assert_equal(qp.re_qpe + 1j * qp.imag_qpe, qp.qpe)
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
        qplist = ncfile.reader.read_qplist_sk(spin=0, sigma_kpoint=[0, 0, 0], ignore_imag=False)
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

        other_qplist = ncfile.reader.read_qplist_sk(spin=0, sigma_kpoint=ncfile.sigma_kpoints[1])
        qpl_merge = qplist.merge(other_qplist)

        assert all(qp in qpl_merge for qp in qplist)
        assert all(qp in qpl_merge for qp in other_qplist)

        if self.has_matplotlib():
            assert qpl_merge.plot_vs_e0(show=False)

        qplist_spin = ncfile.qplist_spin
        assert len(qplist_spin) == ncfile.nsppol
        if self.has_matplotlib():
            assert ncfile.plot_qps_vs_e0(show=False)

    def test_sigeph_robot(self):
        """Tests for SigEPhRobot."""
        #raise self.SkipTest("Disabled")
        filepaths = [
            abidata.ref_file("diamond_444q_SIGEPH.nc"),
        ]
        with abilab.SigEPhRobot.from_files(filepaths) as robot:
            robot.add_file("same_file", filepaths[0])
            repr(robot); str(robot)
            robot.to_string(verbose=2)
            assert len(robot) == 2

            data = robot.get_dataframe()
            assert "qpe" in data

            # Test plot methods
            if self.has_matplotlib():
                assert robot.plot_selfenergy_conv(spin=0, sigma_kpoint=0, band=0, show=False)
                assert robot.plot_selfenergy_conv(spin=0, sigma_kpoint=0, band=0, sortby="nbsum", hue="nqibz", show=False)
                #assert robot.plot_qp_convergence(show=False)
                #assert robot.plot_qps_vs_e0(show=False)
                try:
                    assert robot.plot_qpgaps_t(show=False)
                    assert robot.plot_qpgaps_t(plot_qpmks=True, show=False)
                except ValueError:
                    # workaround for matplotlib bug
                    pass

                assert robot.plot_qpgaps_convergence(itemp=0, sortby="nbsum", show=False)
                assert robot.plot_qpgaps_convergence(itemp=0, sortby="nbsum", hue="nqibz", show=False)

                assert robot.plot_qpdata_convergence(spin=0, sigma_kpoint=(0, 0, 0), band=3,
                        itemp=-1, show=False)
                assert robot.plot_qpdata_convergence(spin=0, sigma_kpoint=(0, 0, 0), band=3,
                        itemp=0, sortby="nbsum", hue="nqibz", show=False)

                # Test plot_qpfield_vs_e0
                assert robot.plot_qpfield_vs_e0("qpeme0", itemp=1, sortby=None, hue=None,
                        colormap="viridis", show=False)
                assert robot.plot_qpfield_vs_e0("ze0", itemp=1, sortby="nbsum", hue="nqibz",
                        colormap="viridis", show=False)

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))
