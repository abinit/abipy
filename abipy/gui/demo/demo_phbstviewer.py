#!/usr/bin/env python
import abipy.data as abidata

from abipy.gui.wxapps import wxapp_phbstviewer

phbst_filepath = abidata.ref_file("trf2_5.out_PHBST.nc")

wxapp_phbstviewer(phbst_filepath).MainLoop()
