"""Tests for htc.JobFile."""
from __future__ import print_function, division

import warnings

from abipy.core.testing import AbipyFileTest
from abipy.htc.jobfile import JobFile


class TestJobFile(AbipyFileTest):
    """Unit tests for JobFile."""

    def setUp(self):
        self.file = JobFile('MyJob.sh', input='myinput')

    def test_shell_line(self):
        """Check the shell line is printed correctly."""
        self.assertContains("#!/bin/bash")
        self.file.shell = 'csh'
        self.assertContains("#!/bin/csh")

    def test_declaration_lines(self):
        """Check the declaration are printed correctly in bash."""
        self.assertContains("""
        MPIRUN=""
        EXECUTABLE=abinit
        INPUT=myinput
        LOG=log
        STDERR=stderr
        """)

    def test_execution_line(self):
        """Check the execution line is printed correctly in bash."""
        self.assertContains("""
        $MPIRUN $EXECUTABLE < $INPUT > $LOG 2> $STDERR
        """)

    def test_mpirun(self):
        """Check mpirun is set correctly."""
        self.file.mpirun = 'openmpirun -n 2'
        self.assertContains("""
        MPIRUN="openmpirun -n 2"
        """)

    def test_modules(self):
        """Check the modules are loaded correctly."""

        self.file.modules = "single"
        self.assertContains("module load single")

        m1, m2 = 'mod1', 'mod2/version/1.4-b'
        lookfor = """
        module load {0}
        module load {1}
        """.format(m1, m2)
        self.file.modules = m1, m2
        self.assertContains(lookfor)

        self.file.modules = [m1, m2]
        self.assertContains(lookfor)

        self.file.set_modules(m1, m2)
        self.assertContains(lookfor)

        self.file.set_modules([m1, m2])
        self.assertContains(lookfor)

    def test_lines_before(self):
        """Check lines_before are printed correctly."""

        single_line = "A single line."
        self.file.lines_before = single_line
        self.assertContains(single_line)

        l1 = "This is my first line!"
        l2 = "And that is my ${SECOND_LINE}"
        lookfor = '\n'.join([l1, l2])

        self.file.lines_before = l1, l2
        self.assertContains(lookfor)

        self.file.lines_before = [l1, l2]
        self.assertContains(lookfor)

        self.file.set_lines_before(l1, l2)
        self.assertContains(lookfor)

        self.file.set_lines_before([l1, l2])
        self.assertContains(lookfor)

    def test_lines_after(self):
        """Check lines_after are printed correctly."""

        single_line = "A single line."
        self.file.lines_after = single_line
        self.assertContains(single_line)

        l1 = "This is my first line!"
        l2 = "And that is my ${SECOND_LINE}"
        lookfor = '\n'.join([l1, l2])

        self.file.lines_after = l1, l2
        self.assertContains(lookfor)

        self.file.lines_after = [l1, l2]
        self.assertContains(lookfor)

        self.file.set_lines_after(l1, l2)
        self.assertContains(lookfor)

        self.file.set_lines_after([l1, l2])
        self.assertContains(lookfor)

    def test_other_lines(self):
        """Check other_lines are printed correctly."""

        single_line = "A single line."
        self.file.other_lines = single_line
        self.assertContains(single_line)

        l1 = "This is my first line!"
        l2 = "And that is my ${SECOND_LINE}"
        lookfor = '\n'.join([l1, l2])

        self.file.other_lines = l1, l2
        self.assertContains(lookfor)

        self.file.other_lines = [l1, l2]
        self.assertContains(lookfor)

        self.file.set_other_lines(l1, l2)
        self.assertContains(lookfor)

        self.file.set_other_lines([l1, l2])
        self.assertContains(lookfor)


class TestJobFileCSH(AbipyFileTest):
    """Unit tests for JobFile with csh shell."""

    def setUp(self):
        self.file = JobFile('MyJob.sh', input='myinput')
        self.file.shell='csh'

    def test_declaration_lines(self):
        """Check the declaration are printed correctly in csh."""
        self.assertContains("""
        set MPIRUN=""
        set EXECUTABLE=abinit
        set INPUT=myinput
        set LOG=log
        set STDERR=stderr
        """)

    def test_execution_line(self):
        """Check the execution line is printed correctly in csh."""
        self.assertContains("""
        ($MPIRUN $EXECUTABLE < $INPUT > $LOG) >& $STDERR
        """)


if __name__ == "__main__":
    import unittest
    unittest.main()

