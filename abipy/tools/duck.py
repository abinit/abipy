# coding: utf-8
"""Duck-typing tests"""
from __future__ import print_function, division, unicode_literals, absolute_import


def is_intlike(obj):
    """
    True if obj represents an integer (float such as 1.0 are included as well).
    """
    # isinstance(i, numbers.Integral)
    try:
        return int(obj) == obj
    except (ValueError, TypeError):
        return False


#def is_listlike(obj):
#   if isinstance(branch, (list, tuple, np.ndarray)):
