#!/usr/bin/env python
"""
Script to create abinit / abipy tutorials in the current work folder.
"""
from __future__ import unicode_literals, division, print_function

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "May 2014"

import os
import os.path

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

def chose_tutorial():
    """
    display a list of available tutorials, ask for a choise, return a tutorial
    """

class Tutorial(object)
    """
    abstract tutorial object defining interface 
    """
    @abstractproprety
    def description(self):
        """
        short description of the topics covered
        """

    @abstractproperty
    def time(self):
        """
        estimated time to follow this tutorial
        """

    def make(self):
        """
        build the tutorial
        """


class ShortBasicsAndGW(Tutorial):
 


if __name__ == "__main__":
     tutorial = chose_tutorial()
     
     tutorial.make()
