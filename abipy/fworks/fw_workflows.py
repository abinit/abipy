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

    def add_to_db(self):
        lpad = LaunchPad.auto_load()
        lpad.add_wf(self.wf)


class ScfFWWorkflow(AbstractFWWorkflow):
    def __init__(self, structure, pseudos):
        abiinput = ebands_input(structure, pseudos).split_datasets()[0]

        abitask = AbiFireTask(abiinput)

        scf_fw = Firework(abitask)

        self.wf = Workflow([scf_fw])
