# coding: utf-8
"""
Base class for parsers
"""
from __future__ import annotations

import abc
import os

from typing import List


class ParserError(Exception):
    """Exceptions raised by Parser objects."""


class BaseParser(metaclass=abc.ABCMeta):
    """
    Abstract class defining the interface that must be provided
    by the parsers used to extract results from the output file of
    a pseudopotential generator a.k.a. ppgen

    Attributes:

        errors: List of strings with errors reported by the pp generator
        warnings: List of strings with the warnings reported by the pp generator.
    """

    Error = ParserError

    def __init__(self, filepath: str) -> None:
        self.filepath = os.path.abspath(filepath)
        self.run_completed = False
        self._errors = []
        self._warnings = []

    @property
    def errors(self) -> List[str]:
        """
        List of strings with possible errors reported by the generator at run-time.
        """
        return self._errors

    @property
    def warnings(self) -> List[str]:
        """
        List of strings with possible errors reported by the generator at run-time.
        """
        return self._warnings

    @abc.abstractmethod
    def get_results(self):
        """
        Return the most important results in a dictionary.
        """

    @abc.abstractmethod
    def get_input_str(self) -> str:
        """Returns a string with the input file."""
