#!/usr/bin/env python
import os
from abipy.gui.wxapps import wxapp_showfiles

# Show all *.py files located within the current directory
app = wxapp_showfiles(filenames=None, dirpath=os.path.dirname(__file__), walk=True, wildcard="*.py")
app.MainLoop()
