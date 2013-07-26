from __future__ import print_function, division

import collections

from abipy.tools import index

__all__ = [
    "Shells",
]


class Shell(collections.Iterable):

    def __init__(self, items, indices, value):
        self.items = items
        self.indices = indices
        self.value = value

        assert len(self.indices) == len(self.items)

    def __str__(self):
        s = "Shell: value = %f, nitems = %i\n" % (self.value, len(self))
        if prtvol:
            for idx, item in self:
                s += "  idx = %d, item = %s,\n" % (idx, str(item))
        return s

    def __iter__(self):
        return iter(zip(self.indices, self.items))

    def __len__(self):
        return len(self.items)

    def __getitem__(self, slice):
        return self.items[slice]

    def indexitem(self):
        for idx, item in zip(self.indices, self.items):
            yield (idx, item)

#########################################################################################


class Shells(collections.Iterable):
    ATOL_SHELL = 1.0e-6

    def __init__(self, iterable, func=None):

        # Build list of tuple (item, idx)
        unsort = [(item, idx) for idx, item in enumerate(iterable)]

        # Sort unsort, then compute the list of values.
        if func is None:
            key = lambda t: t[0]
        else:
            key = lambda t: func(t[0])

        sort_items = sorted(unsort, key=key)

        sort_values = [key(item) for item in sort_items]

        # Build list of shell instance.
        old_val = sort_values[0]
        _shells = list()
        self.values = [old_val]

        indices, items = [], []
        for (i, val) in enumerate(sort_values):
            item = sort_items[i][0]
            index = sort_items[i][1]
            if (val - old_val) > Shells.ATOL_SHELL:
                self.values.append(val)
                sh =  Shell(items, indices, old_val)
                _shells.append(sh)
                #
                old_val = val
                items = [item]
                indices = [index]
            else:
                indices.append(index)
                items.append(item)

        sh = Shell(items, indices, val)
        _shells.append(sh)

        self._shells = tuple(_shells)

    def __str__(self):
        return "\n".join([str(shell) for shell in self])

    def __len__(self):
        return len(self._shells)

    def __iter__(self):
        return self._shells.__iter__()

    def __getitem__(self, slice):
        return self._shells[slice]

    def get_from_value(self, value, atol=None):
        """Return a shell from its value."""
        if atol is None:
            atol = Shells.ATOL_SHELL

        try:
            idx = index(self.values, value, atol=atol)
            return self._shells[idx]
        except ValueError:
            raise ValueError("%s is not in self.values" % str(value))

#########################################################################################
