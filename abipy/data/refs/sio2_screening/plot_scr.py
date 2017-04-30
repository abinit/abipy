#!/usr/bin/env python
"""
This examples shows how to plot the macroscopic dielectric function
computed in the GW code (optdriver 3)
"""
import abipy.data as abidata
from abipy.abilab import abiopen

with abiopen(abidata.ref_file("sio2_SCR.nc")) as ncfile:
    #print(ncfile)
    # The SCR file contains a structure and electron bands in the IBZ.
    # We can thus use the ebands object to plot bands + DOS.
    #edos = ncfile.ebands.get_edos()
    #ncfile.ebands.plot_with_edos(edos, title="KS energies used to compute the SCR file.")
    #ncfile.plot_emacro(title="Macroscopic dielectric function of $SiO_2$ with local-field effects.")
    #ncfile.plot_eelf(title="Electron Energy Loss Function of $SiO_2$")

    #kpoint = [0.0, 0, 0]
    kpoint = [0.5, 0, 0]
    em1 = ncfile.reader.read_wggmat(kpoint)
    print(em1)
    gvec1 = [1, 0, 0]
    gvec2 = [-1, 0, 0]

    em1.plot_freq(gvec1=gvec1, gvec2=gvec2, waxis="real")
    em1.plot_freq(gvec1=[0, 0, 0], gvec2=[1, 0, 0], waxis="imag")

    em1.plot_gg()
    #em1.plot_gg(cplx_mode="abs", wpos=[0], show=False)
    em1.plot_gg(wpos="imag")
