#!/usr/bin/env python
import abipy.data as abidata

from abipy.gui.pseudos import wxapp_pseudos

dirpath = abidata.pseudo_dir
wxapp_pseudos(dirpath).MainLoop()
