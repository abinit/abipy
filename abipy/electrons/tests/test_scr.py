# coding: utf-8
"""Tests for scr module."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np
import abipy.data as abidata

from abipy.core.gsphere import GSphere
from abipy.core.testing import AbipyTest
from abipy.electrons.scr import _AwggMatrix, ScrFile, InverseDielectricFunction


class AwggMatTest(AbipyTest):

    def test_awggmat_api(self):
        """Testing AwggMat API"""
        ecut = 2
        lattice = np.reshape(np.array([1., 0, 0, 0, 1, 0, 0, 0, 1]), (3, 3))
        kpoint = [0, 0, 0]
        gvecs = np.array([[0, 0, 0], [1, 0, 0]])
        gsphere = GSphere(ecut, lattice, kpoint, gvecs, istwfk=1)
        ng = len(gsphere)
        assert ng == len(gvecs)

        wpoints = [0, 1, 2, 3j, 4j]
        nw = len(wpoints)
        data = np.empty((nw, ng, ng), dtype=np.complex)

        f = _AwggMatrix(wpoints, gsphere, data, inord="C")
        repr(f); str(f)

        assert _AwggMatrix.class_from_netcdf_name("inverse_dielectric_function") is InverseDielectricFunction
        with self.assertRaises(ValueError):
            _AwggMatrix.class_from_netcdf_name("foobar")

        assert f.kpoint == gsphere.kpoint
        assert f.ng == len(gsphere)
        assert f.nw == len(wpoints)
        assert f.nrew == 3 and f.nimw == 2 and f.nw == 5
        assert np.all(f.real_wpoints == [0, 1, 2])
        assert np.all(f.imag_wpoints == [3j, 4j])
        self.assert_equal(f.wggmat_realw, f.wggmat[:f.nrew])
        self.assert_equal(f.wggmat_imagw, f.wggmat[f.nrew:])
        assert f.windex(2) == 2
        assert f.windex(3j) == 3
        assert f.gindex([1, 0, 0]) == 1
        assert f.gindex(0) == 0

        for cplx_mode in ("re", "im", "abs", "angle"):
            assert len(f.latex_label(cplx_mode))

        if self.has_matplotlib():
            assert f.plot_w(gvec1=0, gvec2=None, waxis="real", cplx_mode="re-im", show=False)
            assert f.plot_w(gvec1=[0, 0, 0], gvec2=[1, 0, 0], waxis="imag", cplx_mode="re-im", show=False)


class ScrFileTest(AbipyTest):

    def test_scrfile(self):
        """Testing ScrFile."""
        with ScrFile(abidata.ref_file("sio2_SCR.nc")) as ncfile:
            repr(ncfile); str(ncfile)
            ncfile.to_string(verbose=1)
            assert ncfile.netcdf_name == "inverse_dielectric_function"
            assert ncfile.structure.formula == "Si3 O6"
            assert len(ncfile.kpoints) == 9
            assert ncfile.kpoints[0] == [0, 0, 0]
            self.assert_almost_equal(ncfile.kpoints[8].frac_coords, [1/4, 1/4, 1/3])
            assert len(ncfile.wpoints) == 35
            assert ncfile.nw == 35 and ncfile.nrew == 30 and ncfile.nimw == 5
            # TODO Ha or eV?
            #assert ncfile.ng == 9

            params = ncfile.params
            assert params is ncfile.params
            assert params.inclvkb == 2 and params.nbands_used == 50
            assert params.mbpt_sciss == 0 and params.npwwfn_used == 811

            assert len(ncfile.ebands.kpoints) == 9
            assert ncfile.ebands.kpoints.is_mpmesh
            # FIXME: workaround required because fermie is not filled.
            assert ncfile.ebands.fermie == 0

            #em1 = ncfile.get_em1(kpoint=[0, 0, 0])
            #em_nlf = ncfile.get_emacro_nlf(kpoint=(0, 0, 0))
            #em_lf = ncfile.get_emacro_lf(kpoint=(0, 0, 0))

            kpoint = [0.5, 0, 0]
            em1 = ncfile.reader.read_wggmat(kpoint)
            repr(em1); str(em1)
            assert em1.kpoint == kpoint
            assert em1.nw == 35 and em1.nrew == 30 and em1.nimw == 5
            #assert em1.real_wpoints
            #assert em1.imag_wpoints
            assert em1.windex(em1.real_wpoints[1]) == 1
            assert em1.windex(em1.imag_wpoints[1]) == em1.nrew + 1
            assert em1.gindex([0, 0, -1]) == em1.gindex(2)
            for cplx_mode in ("re", "im", "abs", "angle"):
                str(em1.latex_label(cplx_mode))
            with self.assertRaises(ValueError):
                em1.latex_label("foobar")

            if self.has_matplotlib():
                # ncfile plot methods
                #assert em_nlf.plot(show=False)
                assert ncfile.plot_emacro(show=False)
                assert ncfile.ebands.plot(show=False)
                # em1 plot methods
                assert em1.plot_w(gvec1=0, gvec2=None, waxis="real", cplx_mode="re-im", show=False)
                assert em1.plot_w(gvec1=[0, 0, 0], gvec2=[1, 0, 0], waxis="imag", cplx_mode="re-im", show=False)
                #assert em1.plot_gg(cplx_mode="abs", show=False)
                #assert em1.plot_gg(cplx_mode="re", wpos="all", show=False)

            if self.has_nbformat():
                assert ncfile.write_notebook(nbpath=self.get_tmpname(text=True))

