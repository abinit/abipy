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

from fw_tasks import AbiFireTask, ScfFWTask, RelaxFWTask, NscfFWTask, HybridFWTask
from utility_tasks import FinalCleanUpTask
from fw_utils import SHORT_SINGLE_CORE_SPEC, append_fw_to_wf, get_short_single_core_spec
from abipy.abio.factories import ion_ioncell_relax_input, ebands_input, scf_input
from abipy.abio.factories import HybridOneShotFromGsFactory, ScfFactory
from abipy.abio.inputs import AbinitInput

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

    @staticmethod
    def set_short_single_core_to_spec(spec={}):
        spec = dict(spec)

        qadapter_spec = get_short_single_core_spec()
        spec['mpi_ncpus'] = 1
        spec['_queueadapter'] = qadapter_spec
        return spec

    def add_final_cleanup(self, out_exts=["WFK"]):
        cleanup_fw = Firework(FinalCleanUpTask(out_exts=out_exts), spec=self.set_short_single_core_to_spec(),
                              name=(self.wf.name+"_cleanup")[:15])

        append_fw_to_wf(cleanup_fw, self.wf)

    def add_metadata(self, structure=None, additional_metadata={}):
        metadata = dict(wf_type = self.__class__.__name__)
        if structure:
            composition = structure.composition
            metadata['nsites'] = len(structure)
            metadata['elements'] = [el.symbol for el in composition.elements]
            metadata['reduced_formula'] = composition.reduced_formula

        metadata.update(additional_metadata)

        self.wf.metadata.update(metadata)

    def get_reduced_formula(self, input):
        structure = None
        try:
            if isinstance(input, AbinitInput):
                structure = input.structure
            elif 'structure' in input.kwargs:
                structure = input.kwargs['structure']
            else:
                structure = input.args[0]
        except Exception as e:
            logger.warning("Couldn't get the structure from the input: {} {}".format(e.__class__.__name__, e.message))

        return structure.composition.reduced_formula if structure else ""


class InputFWWorkflow(AbstractFWWorkflow):
    def __init__(self, abiinput, task_type=AbiFireTask, autoparal=False, spec={}):
        abitask = task_type(abiinput, is_autoparal=autoparal)

        spec = dict(spec)
        if autoparal:
            spec = self.set_short_single_core_to_spec(spec)

        self.fw = Firework(abitask, spec=spec)

        self.wf = Workflow([self.fw])


class ScfFWWorkflow(AbstractFWWorkflow):
    def __init__(self, abiinput, autoparal=False, spec={}):
        abitask = ScfFWTask(abiinput, is_autoparal=autoparal)

        spec = dict(spec)
        if autoparal:
            spec = self.set_short_single_core_to_spec(spec)

        self.scf_fw = Firework(abitask, spec=spec)

        self.wf = Workflow([self.scf_fw])

    @classmethod
    def from_factory(cls, structure, pseudos, kppa=None, ecut=None, pawecutdg=None, nband=None, accuracy="normal",
                     spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                     shift_mode="Monkhorst-Pack", extra_abivars={}, decorators=[], autoparal=False, spec={}):
        abiinput = scf_input(structure, pseudos, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg, nband=nband,
                             accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
                             scf_algorithm=scf_algorithm, shift_mode=shift_mode)
        abiinput.set_vars(extra_abivars)
        for d in decorators:
            d(abiinput)

        return cls(abiinput, autoparal=autoparal, spec=spec)


class RelaxFWWorkflow(AbstractFWWorkflow):
    def __init__(self, ion_input, ioncell_input, autoparal=False, spec={}):

        spec = dict(spec)
        if autoparal:
            spec = self.set_short_single_core_to_spec(spec)

        ion_task = RelaxFWTask(ion_input, is_autoparal=autoparal)
        self.ion_fw = Firework(ion_task, spec=spec)

        ioncell_task = RelaxFWTask(ioncell_input, is_autoparal=autoparal)
        self.ioncell_fw = Firework(ioncell_task, spec=spec)

        self.wf = Workflow([self.ion_fw, self.ioncell_fw], {self.ion_fw: [self.ioncell_fw]})


class NscfFWWorkflow(AbstractFWWorkflow):
    def __init__(self, scf_input, nscf_input, autoparal=False, spec={}):

        spec = dict(spec)
        if autoparal:
            spec = self.set_short_single_core_to_spec(spec)

        ion_task = ScfFWTask(scf_input, is_autoparal=autoparal)
        self.ion_fw = Firework(ion_task, spec=spec)

        ioncell_task = NscfFWTask(nscf_input, deps={ion_task.task_type: 'DEN'}, is_autoparal=autoparal)
        self.ioncell_fw = Firework(ioncell_task, spec=spec)

        self.wf = Workflow([self.ion_fw, self.ioncell_fw], {self.ion_fw: [self.ioncell_fw]})


class HybridOneShotFWWorkflow(AbstractFWWorkflow):
    def __init__(self, scf_inp, hybrid_input, autoparal=False, spec={}, initialization_info={}):
        rf = self.get_reduced_formula(scf_inp)

        scf_task = ScfFWTask(scf_inp, is_autoparal=autoparal)

        spec = dict(spec)
        spec['initialization_info'] = initialization_info
        if autoparal:
            spec = self.set_short_single_core_to_spec(spec)

        self.scf_fw = Firework(scf_task, spec=spec, name=rf+"_"+scf_task.task_type)

        hybrid_task = HybridFWTask(hybrid_input, is_autoparal=autoparal, deps=["WFK"])

        self.hybrid_fw = Firework(hybrid_task, spec=spec, name=rf+"_"+hybrid_task.task_type)

        self.wf = Workflow([self.scf_fw, self.hybrid_fw], {self.scf_fw: self.hybrid_fw})

    @classmethod
    def from_factory(cls, structure, pseudos, kppa=None, ecut=None, pawecutdg=None, nband=None, accuracy="normal",
                     spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                     shift_mode="Monkhorst-Pack", hybrid_functional="hse06", ecutsigx=None, gw_qprange=1,
                     extra_abivars={}, decorators=[], autoparal=False, spec={}, initialization_info={}):

        scf_fact = ScfFactory(structure=structure, pseudos=pseudos, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg,
                              nband=nband, accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
                              scf_algorithm=scf_algorithm, shift_mode=shift_mode, extra_abivars=extra_abivars,
                              decorators=decorators)

        hybrid_fact = HybridOneShotFromGsFactory(functional=hybrid_functional, ecutsigx=ecutsigx, gw_qprange=gw_qprange,
                                                 decorators=decorators, extra_abivars=extra_abivars)

        return cls(scf_fact, hybrid_fact, autoparal=autoparal, spec=spec, initialization_info=initialization_info)
