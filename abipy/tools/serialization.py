# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Most features of this module has been moved to monty. Please refer to
monty.json and monty.serialization documentation.
"""
from __future__ import annotations

import functools
import json
import pickle
import json

from typing import Any
from pathlib import Path
from monty.json import MontyDecoder, MontyEncoder
from pymatgen.core.periodic_table import Element
from abipy.tools.context_managers import Timer


def pmg_serialize(method):
    """
    Decorator for methods that add MSON serializations keys
    to the dictionary. See documentation of MSON for more details
    """

    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        self = args[0]
        d = method(*args, **kwargs)
        # Add @module and @class
        d["@module"] = type(self).__module__
        d["@class"] = type(self).__name__
        return d

    return wrapper


def json_pretty_dump(obj: Any, filename: str) -> None:
    """
    Serialize obj as a JSON formatted stream to the given filename (
    pretty printing version)
    """
    with open(filename, "w") as fh:
        json.dump(obj, fh, indent=4, sort_keys=4)


class PmgPickler(pickle.Pickler):
    """
    Persistence of External Objects as described in section 12.1.5.1 of
    https://docs.python.org/3/library/pickle.html
    """

    def persistent_id(self, obj: Any):
        """Instead of pickling as a regular class instance, we emit a persistent ID."""
        if isinstance(obj, Element):
            # Here, our persistent ID is simply a tuple, containing a tag and
            # a key
            return type(obj).__name__, obj.symbol
        # If obj does not have a persistent ID, return None. This means obj
        # needs to be pickled as usual.
        return None


class PmgUnpickler(pickle.Unpickler):
    """
    Persistence of External Objects as described in section 12.1.5.1 of
    https://docs.python.org/3/library/pickle.html
    """

    def persistent_load(self, pid):
        """
        This method is invoked whenever a persistent ID is encountered.
        Here, pid is the tuple returned by PmgPickler.
        """
        try:
            type_tag, key_id = pid
        except Exception:
            # Sometimes we get a string such as ('Element', u'C') instead
            # of a real tuple. Use ast to evaluate the expression (much safer
            # than eval).
            import ast
            type_tag, key_id = ast.literal_eval(pid)

        if type_tag == "Element":
            return Element(key_id)

        # Always raises an error if you cannot return the correct object.
        # Otherwise, the unpickler will think None is the object referenced
        # by the persistent ID.
        raise pickle.UnpicklingError(f"unsupported persistent object with pid {pid}")


def pmg_pickle_load(filobj, **kwargs) -> Any:
    """
    Loads a pickle file and deserialize it with PmgUnpickler.

    Args:
        filobj: File-like object
        **kwargs: Any of the keyword arguments supported by PmgUnpickler
    Returns:
        Deserialized object.
    """
    return PmgUnpickler(filobj, **kwargs).load()


def pmg_pickle_dump(obj: Any, filobj, **kwargs):
    """
    Dump an object to a pickle file using PmgPickler.

    Args:
        obj: Object to dump.
        fileobj: File-like object
        **kwargs: Any of the keyword arguments supported by PmgPickler
    """
    return PmgPickler(filobj, **kwargs).dump(obj)


def mjson_load(filepath: str, **kwargs) -> Any:
    """
    Read JSON file in MSONable format with MontyDecoder.
    """
    with open(filepath, "rt") as fh:
        return json.load(fh, cls=MontyDecoder, **kwargs)


def mjson_loads(string: str, **kwargs) -> Any:
    """
    Read JSON string in MSONable format with MontyDecoder.
    """
    return json.loads(string, cls=MontyDecoder, **kwargs)


def mjson_write(obj: Any, filepath: str, **kwargs) -> None:
    """
    Write object to filepath in JSON format using MontyDecoder.
    """
    with open(filepath, "wt") as fh:
        json.dump(obj, fh, cls=MontyEncoder, **kwargs)


class HasPickleIO:
    """
    Mixin class providing pickle IO methods.
    """

    @classmethod
    def pickle_load(cls, workdir, basename=None):
        """
        Reconstruct the object from a pickle file located in workdir.
        """
        filepath = Path(workdir) / f"{cls.__name__}.pickle" if basename is None else Path(workdir) / basename
        with open(filepath, "rb") as fh, Timer(header=f"Reconstructing {cls.__name__} instance from file: {str(filepath)}", footer="") as timer:
            return pickle.load(fh)

    def pickle_dump(self, workdir, basename=None) -> Path:
        """Write pickle file. Return path to file"""
        filepath = Path(workdir) / f"{self.__class__.__name__}.pickle" if basename is None else Path(workdir) / basename
        with open(filepath, "wb") as fh, Timer(header=f"Saving {self.__class__.__name__} instance to file: {str(filepath)}", footer="") as timer:
            pickle.dump(self, fh)
        return filepath
