# coding: utf-8
"""
Firework workflows based on abipy Tasks.
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
from fw_tasks import ScfStrategyFireTask, RelaxStrategyFireTask, NscfStrategyFireTask, AutoparalFireTask, \
    MultiStepRelaxStrategyFireTask, relaxation_methods
from fw_utils import parse_workflow
from pymatgen.io.abinitio.abiobjects import KSampling

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

    def get_task_dirs(self, num_tasks, workdir=None):
        """
        Creates a list of workdirs, one for each task, starting from the original workdir.
        If no workdir is passed it looks in the attributes, and otherwise the running directory is set
        """

        #TODO maybe use a call to map()?
        task_workdirs = []
        for i in range(num_tasks):
            task_workdirs.append(os.path.join(workdir, "t"+str(i)))

        return task_workdirs

    @classmethod
    def create_autoparal_fw(cls, firetask,):
        spec = {'_queueadapter': {'ntasks': 1, 'walltime': '00:10:00'}}
        autoparal_task = AutoparalFireTask(firetask)
        autoparal_fw = Firework(autoparal_task, spec=spec)
        return autoparal_fw


class ScfFWWorkflow(AbstractFWWorkflow):
    """
    Workflow for atomic position relaxations.
    The unit cell parameters are kept fixed
    """

    def __init__(self, structure, pseudos, ksampling=1000, accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 use_symmetries=True, autoparal=False, **extra_abivars):
        """
        Args:

        """

        self.task = ScfStrategyFireTask(structure=structure, pseudos=pseudos, ksampling=ksampling,
                                        accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
                                        scf_algorithm=scf_algorithm, use_symmetries=use_symmetries, **extra_abivars)
        fw = Firework(self.task)

        if autoparal:
            autoparal_fw = self.create_autoparal_fw(self.task)
            self.wf = Workflow([autoparal_fw, fw], {autoparal_fw: [fw]})
        else:
            self.wf = Workflow([fw])


class RelaxAtomsFWWorkflow(AbstractFWWorkflow):
    """
    Workflow for atomic position relaxations.
    The unit cell parameters are kept fixed
    """

    def __init__(self, structure, pseudos, ksampling=1000, accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 autoparal=False, **extra_abivars):
        """
        Args:

        """

        self.task = RelaxStrategyFireTask(structure=structure, pseudos=pseudos, ksampling=ksampling,
                                          relax_algo="atoms_only", accuracy=accuracy, spin_mode=spin_mode,
                                          smearing=smearing, charge=charge, scf_algorithm=scf_algorithm,
                                          **extra_abivars)
        fw = Firework(self.task)

        if autoparal:
            autoparal_fw = self.create_autoparal_fw(self.task)
            self.wf = Workflow([autoparal_fw, fw], {autoparal_fw: [fw]})
        else:
            self.wf = Workflow([fw])

    
class RelaxFWWorkflow(AbstractFWWorkflow):
    """
    Workflow for structural relaxations. The first task relaxes the atomic position
    while keeping the unit cell parameters fixed. The second task uses the final
    structure to perform a structural relaxation in which both the atomic positions
    and the lattice parameters are optimized.
    """
    def __init__(self, structure, pseudos, ksampling=1000, accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 autoparal=False, **extra_abivars):
        """
        Args:

        """

        # Create the ionic relaxation fw
        self.ion_task = RelaxStrategyFireTask(
            structure=structure, pseudos=pseudos, ksampling=ksampling, relax_algo="atoms_only", accuracy=accuracy,
            spin_mode=spin_mode, smearing=smearing, charge=charge, scf_algorithm=scf_algorithm, **extra_abivars)

        ion_fw = Firework(self.ion_task)

        # Create the ionic+cell relaxation fw
        self.ioncell_task = RelaxStrategyFireTask(
            structure=structure, pseudos=pseudos, ksampling=ksampling, relax_algo="atoms_and_cell", accuracy=accuracy,
            spin_mode=spin_mode, smearing=smearing, charge=charge, scf_algorithm=scf_algorithm,
            deps={id(self.ion_task): "WFK structure"}, **extra_abivars)

        ioncell_fw = Firework(self.ioncell_task)

        # Create the workflow
        if autoparal:
            autoparal_fw = self.create_autoparal_fw(self.ion_task)
            self.wf = Workflow([autoparal_fw, ion_fw, ioncell_fw],
                               {autoparal_fw: [ion_fw], ion_fw: [ioncell_fw]})
        else:
            self.wf = Workflow([ion_fw, ioncell_fw], {ion_fw: [ioncell_fw]})


class MultiStepRelaxFWWorkflow(AbstractFWWorkflow):
    """
    Workflow for structural relaxations performed in a multi-step procedure.
    Using a small number of maximum relaxation steps the overall process is
    automatically split in as many FWs are needed.
    """
    def __init__(self, structure, pseudos, ksampling=1000, relax_algo="atoms_only", accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 autoparal=False, max_restart=10, **extra_abivars):
        """
        Args:

        """

        task = MultiStepRelaxStrategyFireTask(structure=structure, pseudos=pseudos, ksampling=ksampling,
                                              relax_algo=relax_algo, accuracy=accuracy, spin_mode=spin_mode,
                                              smearing=smearing, charge=charge, scf_algorithm=scf_algorithm,
                                              deps={}, additional_steps=max_restart, **extra_abivars)

        fw = Firework(task)

        # Create the workflow
        if autoparal:
            autoparal_fw = self.create_autoparal_fw(task)
            self.wf = Workflow([autoparal_fw, fw], {autoparal_fw: [fw]})
        else:
            self.wf = Workflow([fw])


class BandStructureFWWorkflow(AbstractFWWorkflow):
    """
    Workflow to calculate the band structure.
    Can perform a relaxation or just a scf depending on the relax_algo parameter
    """
    def __init__(self, structure, pseudos, ksampling=1000, relax_algo="atoms_and_cell", accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 use_symmetries=True, ksampling_band=8, nscf_bands=None, nscf_algorithm=None, autoparal=False,
                 **extra_abivars):

        # workflow to obtain the density
        if relax_algo in relaxation_methods.keys():
            if relax_algo == "atoms_and_cell":
                dens_wf = RelaxFWWorkflow(structure=structure, pseudos=pseudos, ksampling=ksampling, accuracy=accuracy,
                                          spin_mode=spin_mode, smearing=smearing, charge=charge,
                                          scf_algorithm=scf_algorithm, autoparal=autoparal, **extra_abivars)
                last_task = dens_wf.ioncell_task
            else:
                dens_wf = RelaxAtomsFWWorkflow(structure=structure, pseudos=pseudos, ksampling=ksampling,
                                               accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
                                               scf_algorithm=scf_algorithm, autoparal=autoparal, **extra_abivars)
                last_task = dens_wf.task
            deps = {id(last_task): 'DEN structure'}
        else:
            if relax_algo not in (None, 'scf', 'scf_only', 'no_relax'):
                logger.warning('relax_algo value "%s" not recognized. Perorming a scf calculation.' % relax_algo)

            dens_wf = ScfFWWorkflow(structure=structure, pseudos=pseudos, ksampling=ksampling, accuracy=accuracy,
                                    spin_mode=spin_mode, smearing=smearing, charge=charge, scf_algorithm=scf_algorithm,
                                    use_symmetries=use_symmetries, autoparal=autoparal, **extra_abivars)
            last_task = dens_wf.task
            deps = {id(last_task): 'DEN'}

        # Create the ksampling based on high symmetry lines of the structure
        if isinstance(ksampling_band, int):
            # ksampling_band = KSampling.automatic_density(structure, 100)
            ksampling_band = KSampling.path_from_structure(ksampling_band, structure)

        # Create the nscf FW
        nscf_task = NscfStrategyFireTask(scf_task=last_task, ksampling=ksampling_band, nscf_bands=nscf_bands,
                                         nscf_algorithm=nscf_algorithm, deps=deps, **extra_abivars)
        nscf_fw = Firework(nscf_task)

        # Expand the first part to get the proper list of FWs and links
        fws, links = parse_workflow([dens_wf.wf, nscf_fw], {dens_wf.wf: [nscf_fw]})

        self.wf = Workflow(fws, links)

