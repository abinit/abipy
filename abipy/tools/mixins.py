"""Collection of mixin classes."""
from __future__ import print_function, division

import collections


class MappingMixin(collections.Mapping):
    """
    Mixin class implementing the mapping protocol. Useful to avoid boilerplate code if you want
    to define a object that behaves as a Mapping but without inheriting from dict or similar classes
    because you don't want to expose/support all the methods of dict.

    Client code must initialize a Mapping object either in new or in init and bound it to _mapping_mixin_
    The implemention of the Mapping interface is delegated to _mapping_mixin_

    .. Example:

    >>> class Foo(MappingMixin):
    ...     def __init__(self, attr, **kwargs):
    ...         self._mapping_mixin_ = kwargs
    ...         self.attr = attr
    >>> obj = Foo(attr=1, spam=2)
    >>> obj.attr, obj["spam"]
    (1, 2)
    >>> obj.pop("spam")
    2
    >>> len(obj), "spam" in obj
    (0, False)
    """
    def __len__(self):
        return len(self._mapping_mixin_)

    def __iter__(self):
        return self._mapping_mixin_.__iter__()

    def __getitem__(self, key):
        return self._mapping_mixin_[key]

    def __setitem__(self, key, value):
        self._mapping_mixin_[key] = value

    def __contains__(self, key):
        return key in self._mapping_mixin_

    def keys(self):
        return self._mapping_mixin_.keys()

    def items(self):
        return self._mapping_mixin_.items()

    def get(self, value, default=None):
        try:
            return self[value]
        except KeyError:
            return default

    def pop(self, k, *d):
        """
        D.pop(k[,d]) -> v, remove specified key and return the corresponding value.
        If key is not found, d is returned if given, otherwise KeyError is raised
        """
        if d:
            return self._mapping_mixin_.pop(k, d[0])
        else:
            return self._mapping_mixin_.pop(k)

    def update(self, e, **f):
        """
        D.update([E, ]**F) -> None.  Update D from dict/iterable E and F.
        If E present and has a .keys() method, does:     for k in E: D[k] = E[k]
        If E present and lacks .keys() method, does:     for (k, v) in E: D[k] = v
        In either case, this is followed by: for k in F: D[k] = F[k]
        """
        self._mapping_mixin_.update(e, **f)


if __name__ == "__main__":
    import unittest
    unittest.main()
