#!/usr/bin/env python

import os
import abipy
import abipy.gui.apps as apps

filenames = [os.path.join(abipy.get_datadir(), f) for f in  ["t01.out", "t02.out"]]
apps.awx_events(filenames)
