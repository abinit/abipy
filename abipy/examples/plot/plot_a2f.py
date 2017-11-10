#!/usr/bin/env python
r"""
Eliashberg function
===================

This example shows how to plot the Eliashberg function a2F(w)
and the e-ph coupling strenght in metals.
"""
from abipy import abilab
import abipy.data as abidata

phdos_path = abidata.ref_file("al_161616q_PHDOS.nc")

with abilab.abiopen(abidata.ref_file("al_888k_161616q_EPH.nc")) as eph:
    print(eph)
    #eph.phbands.plot()

    #eph.a2f_qintp.plot()
    #with_lambda = False
    #fig = eph.a2f_qcoarse.plot_nuterms(with_lambda=with_lambda, show=False)
    #eph.a2f_qintp.plot_nuterms(axmat=fig.axes, with_lambda=with_lambda)

    #eph.plot()
    #eph.plot_eph_strength(what="lambda")
    eph.plot_with_a2f(phdos=phdos_path)

    # Grid with 3 plots (a2F, F, a2F) with F taken from PHDOS file.
    #eph.a2f_qintp.plot_a2(phdos_path)
