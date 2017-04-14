#!/usr/bin/env python
"""Tests for visualizer module"""
from __future__ import print_function, division, unicode_literals, absolute_import

#import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.iotools.visualizer import Visualizer, Xcrysden, Vesta, V_Sim


class TestVisualizer(AbipyTest):

    def test_visualizers(self):
        print("Available visualizers:")
        for vclass in Visualizer.get_available():
            print(vclass)

        assert Xcrysden.support_ext("xsf") and Xcrysden.support_ext(".xsf")

        assert Vesta is Visualizer.from_name("vesta")
        assert Vesta.support_ext("xsf") and "xsf" in Vesta.supported_extensions()
        #assert Vesta.from_file("foo.xsf")

        assert V_Sim.support_ext("xsf")

        for cls in [Xcrysden, Vesta, V_Sim]:
            visu = cls("foo.xsf")
            assert callable(visu)
            str(visu)
            repr(visu)
            # cmdarg is a string?
            assert visu.cmdarg + " "
