#!/usr/bin/env python
import abipy.data as abidata
from abipy.gui.wxapps import wxapp_mdfviewer

mdf_filenames = abidata.ref_file("tbs_4o_DS2_MDF.nc")

wxapp_mdfviewer(mdf_filenames).MainLoop()
