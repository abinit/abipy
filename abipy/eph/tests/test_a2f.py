"""Tests for a2f module."""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest


class A2fFileTest(AbipyTest):

    def test_a2fnc_file(self):
        """Tests for A2fFile."""
        ncfile = abilab.abiopen(abidata.ref_file("al_888k_161616q_A2F.nc"))
        repr(ncfile); str(ncfile)
        assert ncfile.to_string(verbose=2)
        assert ncfile.params["nspinor"] == ncfile.nspinor
        assert "eph_fsewin" in ncfile.params
        assert ncfile.structure.formula == "Al1" and len(ncfile.structure) == 1
        # Ebands
        assert ncfile.nsppol == 1 and ncfile.nspden == 1 and ncfile.nspinor == 1
        assert ncfile.ebands.kpoints.is_ibz
        #self.assert_equal(ncfile.ebands.kpoints.ksampling.mpdivs, [8, 8, 8])
        self.assert_equal(ncfile.ebands.kpoints.ksampling.mpdivs, [12, 12, 12])
        # Phbands
        assert ncfile.phbands.qpoints.is_path
        assert ncfile.phbands.qpoints.ksampling is None
        assert ncfile.phbands.has_linewidths

        # Test ksampling
        ksamp = ncfile.reader.read_ksampling_info()
        assert ksamp.is_mesh and not ksamp.is_path
        assert ksamp.has_diagonal_kptrlatt
        self.assert_equal(ksamp.mpdivs, [12, 12, 12])
        self.assert_equal(ksamp.shifts.ravel(), [0, 0, 0])
        assert ksamp.to_string(title="Ksampling")

        # Test edos
        # TODO
        #ncfile.edos
        #if self.has_matplotlib():
            #assert ncfile.edos.plot(show=False)

        # Test A2f(w) function.
        a2f = ncfile.a2f_qcoarse
        assert ncfile.get_a2f_qsamp("qcoarse") is a2f
        repr(a2f); str(a2f)
        a2f.to_string(verbose=2)
        assert a2f.nsppol == ncfile.nsppol and a2f.nmodes == 3 * len(ncfile.structure)
        assert a2f.iw0 == 0
        #assert a2f.mesh
        #assert a2f.values_spin
        m1 = a2f.get_moment(n=1)
        #self.assert_almost_equal(m1/2, a2f.get_moment(n=1, spin=0))
        #self.assert_almost_equal(self.lambda_iso, )
        #self.assert_almost_equal(self.omega_log, )
        tc = a2f.get_mcmillan_tc(mustar=0.1)
        #self.assert_almost_equal(tc, )
        mustar = a2f.get_mustar_from_tc(tc)
        self.assert_almost_equal(mustar, 0.1)
        #self.assert_almost_equal(a2f.get_mcmillan_tc(mustar), tc)

        assert not ncfile.has_a2ftr
        assert ncfile.a2ftr_qcoarse is None
        assert ncfile.a2ftr_qintp is None
        assert ncfile.get_a2ftr_qsamp("qcoarse") is ncfile.a2ftr_qcoarse
        phdos_path = abidata.ref_file("al_161616q_PHDOS.nc")

        if self.has_matplotlib():
            # Test A2f plot methods
            assert ncfile.a2f_qcoarse.plot(exchange_xy=True, fontsize=8, units="meV", show=False)
            assert ncfile.a2f_qcoarse.plot_with_lambda(show=False)
            assert ncfile.a2f_qcoarse.plot_tc_vs_mustar(show=False)
            assert ncfile.a2f_qintp.plot_a2(phdos_path, show=False)
            assert ncfile.a2f_qintp.plot_nuterms(show=False)
            assert ncfile.plot_a2f_interpol(show=False)

            # Test A2fFile plot methods.
            assert ncfile.plot(show=False)
            assert ncfile.plot_with_a2f(show=False)
            assert ncfile.plot_eph_strength(show=False)
            assert ncfile.plot_with_a2f(qsamp="qcoarse", show=False)
            assert ncfile.plot_with_a2f(phdos=phdos_path, show=False)

        if self.has_nbformat():
            ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

        ncfile.close()


class A2fRobotTest(AbipyTest):

    def test_a2f_robot(self):
        """Test A2fRobot."""
        files = abidata.ref_files(
                "al_888k_161616q_A2F.nc",
                #"al_888k_161616q_A2F.nc",
        )
        with abilab.A2fRobot.from_files(files[0]) as robot:
            robot.add_file("same_a2f", files[0])
            assert len(robot) == 2
            repr(robot); str(robot)
            robot.to_string(verbose=2)
            #assert [t[2] for t in robot.sortby("nkpt")] == [10, 60, 182]

            df_params = robot.get_params_dataframe()
            self.assert_equal(df_params["nspden"].values, 1)

            data = robot.get_dataframe(with_geo=True)
            assert "lambda_qcoarse" in data and "omegalog_qintp" in data

            # Mixin
            phbands_plotter = robot.get_phbands_plotter()
            data = robot.get_phbands_dataframe()
            assert "min_freq" in data
            assert np.all(data["nmodes"].values == 3)

            # Test plot methods
            if self.has_matplotlib():
                assert phbands_plotter.boxplot(show=False)

                # Test wrappers provided by RobotWithPhbands
                assert robot.combiplot_phbands(show=False)
                assert robot.gridplot_phbands(show=False)
                assert robot.boxplot_phbands(show=False)
                assert robot.combiboxplot_phbands(show=False)

                # Test A2fRobot plot methods
                assert robot.plot_lambda_convergence(sortby=None, hue=None, show=False)
                assert robot.plot_lambda_convergence(what="gamma", sortby=None, hue="nkpt", show=False)
                assert robot.plot_a2f_convergence(show=False)
                assert robot.plot_a2f_convergence(hue="nkpt", show=False)
                assert robot.plot_a2fdata_convergence(show=False, sortby=None, hue="nkpt")
                assert robot.gridplot_a2f(show=False)

                #assert robot.plot_a2ftr_convergence(show=False, sortby=None, hue="nkpt")

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))
