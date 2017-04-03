# coding: utf-8
"""Test for output files"""
from __future__ import unicode_literals, division, print_function

import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.abio.outputs import AbinitOutputFile


class AbinitOutputTest(AbipyTest):

    def test_gs_output(self):
        """Testing AbinitOutputFile with GS calculation."""
        abo_path = abidata.ref_file("refs/si_ebands/run.abo")
        with AbinitOutputFile(abo_path) as gs_abo:
            print(gs_abo)
            print(gs_abo.events)
            gs_cycle = gs_abo.next_gs_scf_cycle()
            assert gs_cycle is not None
            if self.has_matplotlib():
                gs_cycle.plot(show=False)
            gs_abo.seek(0)
            assert gs_abo.next_d2de_scf_cycle() is None

            timer = gs_abo.get_timer()
            assert len(timer) == 1
            str(timer.summarize())

            if self.has_matplotlib():
                gs_abo.compare_gs_scf_cycles([abo_path], show=False)
                timer.plot_all()

            if self.has_nbformat():
                gs_abo.write_notebook(nbpath=self.get_tmpname(text=True))
                timer.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_ph_output(self):
        """Testing AbinitOutputFile with phonon calculations."""
        # TODO: File is missing
        #abo_path = abidata.ref_file("refs/si_ebands.run.abo"))
        #with AbinitOutputFile(abo_path) as gs_abo:
        #    print(ph_abo)
        #    print(gs_abo.events)
        #    assert ph_abo.next_gs_scf_cycle() is None

        #    ph_abo.seek(0)
        #    ph_cycle = ph_abo.next_d2de_scf_cycle()
        #    assert ph_cycle is not None
        #    if self.has_matplotlib():
        #        ph_cycle.plot(show=False)
        #        gs_abo.compare_d2de_cycles([abo_path], show=False)

        #    if self.has_nbformat():
        #        ph_abo.write_notebook(nbpath=self.get_tmpname(text=True))
