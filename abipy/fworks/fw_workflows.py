# coding: utf-8
"""
Firework workflows
"""
from __future__ import print_function, division, unicode_literals

try:
    from fireworks.core.launchpad import LaunchPad
    from fireworks.core.firework import Firework, Workflow
except ImportError:
    LaunchPad, Workflow, Firework = 3 * [object]

import abc
import six
import os
import logging
import sys
from fw_tasks import AbiFireTask
from abipy.htc.factories import ion_ioncell_relax_input, ebands_input

# logging.basicConfig()
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stdout))

@six.add_metaclass(abc.ABCMeta)
class AbstractFWWorkflow():
    """
    Abstract Workflow class.
    """

    def add_to_db(self, lpad=None):
        if not lpad:
            lpad = LaunchPad.auto_load()
        return lpad.add_wf(self.wf)


class InputFWWorkflow(AbstractFWWorkflow):
    def __init__(self, abiinput):
        abitask = AbiFireTask(abiinput)

        self.fw = Firework(abitask)

        self.wf = Workflow([self.fw])


class ScfFWWorkflow(AbstractFWWorkflow):
    def __init__(self, structure, pseudos):
        abiinput = ebands_input(structure, pseudos).split_datasets()[0]

        abitask = AbiFireTask(abiinput)

        self.scf_fw = Firework(abitask)

        self.wf = Workflow([self.scf_fw])
