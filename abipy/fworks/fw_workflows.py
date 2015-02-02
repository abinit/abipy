# coding: utf-8
"""
Firework workflows based on abipy Tasks.
"""
from __future__ import print_function, division, unicode_literals

import abc
import six
import os

from fireworks.core.launchpad import LaunchPad
from fireworks.core.firework import Firework, Workflow
from fw_tasks import * #RelaxFWTask, MultiStepRelaxFWTask, AutoparalFWTask, AbiFireTask
from pymatgen.io.abinitio.tasks import Dependency, RelaxTask


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
    def create_autoparal_fw(self, firetask, max_ncpus=None):
        spec = {'_queueadapter': {'ntasks': 1, 'walltime': '00:10:00'}}
        autoparal_task = AutoparalFireTask(firetask, max_ncpus)
        autoparal_fw = Firework(autoparal_task, spec=spec)
        return autoparal_fw



class RelaxAtomsFWWorkflow(AbstractFWWorkflow):
    """
    Workflow for atomic position relaxations.
    The unit cell parameters are kept fixed
    """

    #TODO check if manager is useful at this point
    def __init__(self, structure, pseudos, ksampling=1000, accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 autoparal=False, max_ncpus=32, **extra_abivars):
        """
        Args:

        """

        task = RelaxStrategyFireTask(structure=structure, pseudos=pseudos, ksampling=ksampling, relax_algo="atoms_only",
                                     accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
                                     scf_algorithm=scf_algorithm, **extra_abivars)
        fw = Firework(task)

        if autoparal:
            autoparal_fw = self.create_autoparal_fw(task, max_ncpus)
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
                 autoparal=False, max_ncpus=32, **extra_abivars):
        """
        Args:

        """

        # Create the ionic relaxation fw
        ion_task = RelaxStrategyFireTask(structure=structure, pseudos=pseudos, ksampling=ksampling,
                                         relax_algo="atoms_only", accuracy=accuracy, spin_mode=spin_mode,
                                         smearing=smearing, charge=charge,scf_algorithm=scf_algorithm, **extra_abivars)

        ion_fw = Firework(ion_task)

        # Create the ionic+cell relaxation fw
        ioncell_task = RelaxStrategyFireTask(structure=structure, pseudos=pseudos, ksampling=ksampling,
                                             relax_algo="atoms_and_cell", accuracy=accuracy, spin_mode=spin_mode,
                                             smearing=smearing, charge=charge,scf_algorithm=scf_algorithm,
                                             deps={id(ion_task): "WFK structure"},**extra_abivars)

        ioncell_fw = Firework(ioncell_task)

        # Create the workflow
        if autoparal:
            autoparal_fw = self.create_autoparal_fw(ion_task, max_ncpus)
            self.wf = Workflow([autoparal_fw, ion_fw, ioncell_fw],
                               {autoparal_fw: [ion_fw], ion_fw: [ioncell_fw]})
        else:
            self.wf = Workflow([ion_fw, ioncell_fw], {ion_fw: [ioncell_fw]})


#TODO fix after refactoring
class MultiStepRelaxFWWorkflow(AbstractFWWorkflow):
    """
    Workflow for structural relaxations performed in a multi-step procedure.
    Using a small number of maximum relaxation steps the overall process is
    automatically split in as many FWs are needed.
    """
    #TODO check if manager is useful at this point
    def __init__(self, input, workdir=None, manager=None):
        """
        Args:
            ion_input:
                Input for the relaxation of the ions (cell is fixed)
            ioncell_input:
                Input for the relaxation of the ions and the unit cell.
            workdir:
                Working directory.
            manager:
                `TaskManager` object.
        """

        task = MultiStepRelaxFWTask.from_input(input, workdir, manager)

        fw = Firework(task)

        self.wf = Workflow([fw])