# coding: utf-8
"""Test for output files"""
from __future__ import unicode_literals, division, print_function

import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.abio.outputs import AbinitOutputFile, AbinitLogFile


class AbinitLogFileTest(AbipyTest):

    def test_abinit_logfile(self):
        """"Testing AbinitLogFile."""
        log_path = abidata.ref_file("refs/abinit.log")
        with AbinitLogFile(log_path) as abilog:
            str(abilog)
            assert len(abilog.events) == 2
            if self.has_nbformat():
                abilog.write_notebook(nbpath=self.get_tmpname(text=True))


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
        abo_path = abidata.ref_file("refs/gs_dfpt.abo")
        with AbinitOutputFile(abo_path) as abo:
             str(abo)
             gs_cycle = abo.next_gs_scf_cycle()
             assert gs_cycle is not None
             #abo.seek(0)
             ph_cycle = abo.next_d2de_scf_cycle()
             assert ph_cycle is not None
             if self.has_matplotlib():
                ph_cycle.plot(show=False)
                abo.compare_d2de_scf_cycles([abo_path], show=False)

             if self.has_nbformat():
                abo.write_notebook(nbpath=self.get_tmpname(text=True))
