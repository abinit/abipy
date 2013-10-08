"""Tests for htc.InputFile."""
import warnings

from abipy.htc import *
from abipy.tests import *

# =========================================================================== #

class TestAbinitVariable(AbipyFileTest):
    """Unit tests for AbinitVariable."""

    def setUp(self):
        self.file = AbinitVariable('ecut', 10.0)

    def test_varnames(self):
        """Test printing of variables name."""
        self.file.name = 'ecut1'
        self.assertContains('ecut1')

        self.file.name = 'ecut_s'
        self.assertContains('ecut:')

        self.file.name = 'ecut_i'
        self.assertContains('ecut+')

        self.file.name = 'ecut_a'
        self.assertContains('ecut?')

        self.file.name = 'ecut_s2'
        self.assertContains('ecut:2')

        self.file.name = 'ecut3_a'
        self.assertContains('ecut3?')

        self.file.name = 'ecut_s_a'
        self.assertContains('ecut:?')

    def test_scalar_values(self):
        """Test printing of scalar variables."""
        self.file.value = 11.5
        self.assertContains('ecut 11.5')

        self.file.value = 10
        self.assertContains('ecut 10')

        self.file.value = '*1'
        self.assertContains('ecut *1')

        self.file.value = None
        self.assertEmpty()

        self.file.value = ''
        self.assertEmpty()


class TestInputFile(AbipyFileTest):
    """Unit tests for InputFile."""

    def setUp(self):
        self.file = InputFile('MyDir/MyName.in')

    def test_set_variable_attribute(self):
        """Test setting variables by attribute."""
        self.file.ecut = 10.
        self.assertContains('ecut 10.')

    def test_set_variable_function(self):
        """Test setting variables with set_variable."""
        self.file.set_variable('ecut', 10.)
        self.assertContains('ecut 10.')

    def test_set_variables_function(self):
        """Test setting variables with set_variables."""

        self.file.set_variables({'ecut':10., 'nstep':100})
        self.assertContains('nstep 100')
        self.assertContains('ecut 10.')

        self.file.set_variables({'ecut':10., 'nstep':100}, 1)
        self.assertContains('nstep1 100')
        self.assertContains('ecut1 10.')

