#!/usr/bin/env python

import os
import abipy

import abipy.gui.wxapps as wxapps

wxapps.wxapp_dirbrowser(abipy.get_datadir()).MainLoop()

