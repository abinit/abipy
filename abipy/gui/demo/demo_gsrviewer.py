#!/usr/bin/env python
import abipy.data as abidata

from abipy.gui.wxapps import wxapp_gsrviewer

gsr_filename = abidata.ref_file("si_nscf_GSR.nc")

wxapp_gsrviewer(gsr_filename).MainLoop()
