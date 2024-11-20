# coding: utf-8
"""Tests for duck module."""
import numpy as np

from abipy.core.testing import AbipyTest
from abipy.tools import duck


class DuckTest(AbipyTest):

    def test_is_string(self):
        """Testing is_string."""
        assert duck.is_string("hello")
        assert not duck.is_string(1)
        assert not duck.is_string([1, 2, 3])

    def test_is_number_like(self):
        """Testing is_number_like."""
        from numbers import Number
        from decimal import Decimal
        from fractions import Fraction
        is_number_like = duck.is_number_like
        assert is_number_like(2)
        assert is_number_like(2.0)
        assert is_number_like(Decimal('2.0'))
        assert is_number_like(complex(2,0))
        assert is_number_like(Fraction(2,1))
        assert not is_number_like('2')
        assert not is_number_like({})

    def test_is_intlike(self):
        """Testing is_intlike."""
        assert duck.is_intlike(1)
        assert duck.is_intlike(1.0)
        assert not duck.is_intlike(1.3)
        assert not duck.is_intlike("1")
        assert not duck.is_intlike([1, 2, 3])

        for t in (int, np.int64, np.float64):
            assert duck.is_intlike(t(123))

        assert duck.is_intlike(complex(123))
        assert not duck.is_intlike(complex(123.2))
        assert not duck.is_intlike(123 + 1j * 2)

    def test_is_listlike(self):
        """Testing is_listlike."""
        #if isinstance(branch, (list, tuple, np.ndarray)):
        assert duck.is_listlike([1, 2, 3])
        assert duck.is_listlike((1,))
        assert not duck.is_listlike({1, 2, 3})
        assert duck.is_listlike(np.zeros(3))
        d = {"hello": 1, "world": 2}
        #assert not duck.is_listlike(d)
        #assert duck.is_listlike(range(1, 3))

    def test_list_ints(self):
        """Testing list_ints."""
        assert duck.list_ints(1) == [1]
        assert duck.list_ints([1]) == [1]
        assert duck.list_ints([1, 2, 3]) == [1, 2, 3]

    def test_getattrd(self):
        """Testing getattrd and hasattrd."""
        getattrd, hasattrd = duck.getattrd, duck.hasattrd
        with self.assertRaises(AttributeError): getattrd(int, 'a')
        assert getattrd(int, 'a', default=None) is None
        assert getattrd(int, '__class__.__name__') == "type"
        assert hasattrd(int, '__class__.__name__')
        assert not hasattrd(int, 'foobar.__name__')
