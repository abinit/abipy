"""Tests for lessons"""
from __future__ import print_function, division

import os

from abipy.core.testing import *


class TestLessons(AbipyTest):
    """Unit tests for lessons."""

    @staticmethod
    def assert_lesson_object(lesson):
        """Helper function to test `Lesson` protocol."""
        assert len(str(lesson))
        assert lesson.abipy_string + " "
        assert lesson.comline_string + " "
        assert os.path.exists(lesson.pyfile)

    def test_lesson_base1(self):
        """Testing lesson_base1."""
        from abipy.lessons.lesson_base1 import build_flow, analyze_flow
        #self.assert_lesson_object(Lesson())
        flow = build_flow()
        flow.make_scheduler().start()
        analyze_flow(flow)
        flow.rmtree()

    def test_lesson_bse(self):
        """Testing lesson_bse."""
        from abipy.lessons.lesson_bse import make_scf_nscf_bse_inputs

        scf_input, nscf_input, bse_input = make_scf_nscf_bse_inputs(
            ngkpt=(6, 6, 6), ecut=6, ecuteps=3, mdf_epsinf=12.0, mbpt_sciss="0.8 eV")

        bse_input.abivalidate()

    def test_lesson_dfpt(self):
        """Testing lesson_dfpt."""
        from abipy.lessons.lesson_dfpt import make_scf_input
        scf_input = make_scf_input(ecut=2, ngkpt=(4, 4, 4))
        scf_input.abivalidate()

    def test_lesson_dos_bands(self):
        """Testing lesson_dos_bands."""
        from abipy.lessons.lesson_dos_bands import Lesson
        lesson = Lesson()
        self.assert_lesson_object(lesson)

        flow = lesson.make_flow()
        flow.make_scheduler().start()
        #flow.build_and_pickle_dump()
        #lesson.setup()
        #lesson.analyze(flow)
        flow.rmtree()

    def test_lesson_ecut_convergence(self):
        """Testing lesson_ecut_convergence."""
        from abipy.lessons.lesson_ecut_convergence import Lesson
        lesson = Lesson()
        self.assert_lesson_object(lesson)

        flow = lesson.make_ecut_flow()
        flow.make_scheduler().start()
        #flow.build_and_pickle_dump()
        #lesson.setup()
        #lesson.analyze(flow)
        flow.rmtree()

    def test_lesson_g0w0(self):
        """Testing lesson_g0w0."""
        from abipy.lessons.lesson_g0w0 import Lesson
        lesson = Lesson()
        self.assert_lesson_object(lesson)

        flow = lesson.make_flow()
        #flow.build_and_pickle_dump()
        flow.make_scheduler().start()
        #lesson.setup()
        #lesson.analyze(flow)
        flow.rmtree()

    def test_lesson_kpoint_convergence(self):
        """Testing lesson_kpoint_convergence."""
        from abipy.lessons.lesson_kpoint_convergence import Lesson
        lesson = Lesson()
        self.assert_lesson_object(lesson)

        flow = lesson.make_ngkpt_flow()
        flow.make_scheduler().start()
        flow.rmtree()

    def test_lesson_relaxation(self):
        """Testing lesson_relaxation."""
        from abipy.lessons.lesson_relaxation import Lesson
        lesson = Lesson()
        self.assert_lesson_object(lesson)
        flow = Lesson.make_eos_flow()
        flow = Lesson.make_relax_flow()

    #def test_lesson_spin(self):
    #    """Testing lesson_spin."""
    #    from abipy.lessons.lesson_spin import gs_flow
