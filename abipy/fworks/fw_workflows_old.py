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
from fw_tasks_old import NscfStrategyFireTask, AutoparalFireTask, \
    MultiStepRelaxStrategyFireTask, relaxation_methods, AbiFireTask, RelaxFWTask
from fw_utils import parse_workflow
from pymatgen.io.abinitio.abiobjects import KSampling
from pymatgen.io.abinitio.strategies import StrategyWithInput
from abipy.htc.input_factories import ion_ioncell_relax_input
from pymatgen.io.abinitio.tasks import RelaxTask

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
    def create_autoparal_fw(cls, firetask, folder=None):
        spec = {'_queueadapter': {'ntasks': 1, 'walltime': '00:10:00'}}
        if folder:
            spec['_launch_dir'] = folder + '_autoparal'
        autoparal_task = AutoparalFireTask(firetask)
        autoparal_fw = Firework(autoparal_task, spec=spec)
        return autoparal_fw

    autoparal_spec = {'_queueadapter': {'ntasks': 1, 'walltime': '00:10:00'}}



# class ScfFWWorkflow(AbstractFWWorkflow):
#     """
#     Workflow for atomic position relaxations.
#     The unit cell parameters are kept fixed
#     """
#
#     def __init__(self, structure, pseudos, ksampling=1000, accuracy="normal",
#                  spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
#                  use_symmetries=True, autoparal=False, folder=None, **extra_abivars):
#         """
#         Args:
#
#         """
#
#         self.task = ScfStrategyFireTask(structure=structure, pseudos=pseudos, ksampling=ksampling,
#                                         accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
#                                         scf_algorithm=scf_algorithm, use_symmetries=use_symmetries, **extra_abivars)
#
#         spec = {}
#         if folder:
#             spec['_launch_dir'] = os.path.join(folder, 'scf')
#         fw = Firework(self.task, spec=spec)
#
#         if autoparal:
#             autoparal_fw = self.create_autoparal_fw(self.task, spec['_launch_dir'])
#             self.wf = Workflow([autoparal_fw, fw], {autoparal_fw: [fw]})
#         else:
#             self.wf = Workflow([fw])


class RelaxAtomsFWWorkflow(AbstractFWWorkflow):
    """
    Workflow for atomic position relaxations.
    The unit cell parameters are kept fixed
    """

    def __init__(self, structure, pseudos, kppa=1000, ecut=None, pawecutdg=None, accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 autoparal=False, folder=None, **extra_abivars):
        """
        Args:

        """
        abiinput = ion_ioncell_relax_input(structure=structure, pseudos=pseudos, kppa=kppa, ecut=ecut,
                                           pawecutdg=pawecutdg, accuracy=accuracy, spin_mode=spin_mode,
                                           smearing=smearing, charge=charge, scf_algorithm=scf_algorithm)

        abiinput.set_vars(**extra_abivars)

        abitask = RelaxTask(StrategyWithInput(abiinput.split_datasets()[0]))

        self.task = RelaxFWTask(abitask)

        spec = {}
        if folder:
            spec['_launch_dir'] = os.path.join(folder, 'atomic_relax')
        fw = Firework(self.task, spec=spec)

        if autoparal:
            autoparal_fw = self.create_autoparal_fw(self.task, spec['_launch_dir'])
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

    def __init__(self, structure, pseudos, kppa=1000, ecut=None, pawecutdg=None, accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 autoparal=False, folder=None, aim_dilatmx=1.01, **extra_abivars):
        """
        Args:

        """

        abiinput = ion_ioncell_relax_input(structure=structure, pseudos=pseudos, kppa=kppa, ecut=ecut,
                                           pawecutdg=pawecutdg, accuracy=accuracy, spin_mode=spin_mode,
                                           smearing=smearing, charge=charge, scf_algorithm=scf_algorithm)
        abiinput.set_vars(**extra_abivars)
        ion_input, ioncell_input = abiinput.split_datasets()

        abi_ion_task = RelaxTask(StrategyWithInput(ion_input))

        # Create the ionic relaxation fw
        self.ion_task = RelaxFWTask(abi_ion_task, is_autoparal=autoparal)

        spec = {}
        if folder:
            spec['_launch_dir'] = os.path.join(folder, 'atomic_relax')
        if autoparal:
            spec.update(self.autoparal_spec)

        ion_fw = Firework(self.ion_task, spec=spec)

        abi_ioncell_task = RelaxTask(StrategyWithInput(ioncell_input))

        # Create the ionic+cell relaxation fw
        self.ioncell_task = RelaxFWTask(abi_ioncell_task, deps={id(self.ion_task): "structure"}, is_autoparal=autoparal)

        spec = {'create_file': 'final'}
        if folder:
            spec['_launch_dir'] = os.path.join(folder, 'atomic_and_cell_relax')
        if autoparal:
            spec.update(self.autoparal_spec)

        ioncell_fw = Firework(self.ioncell_task, spec=spec)

        # # Create the workflow
        # if autoparal:
        #     autoparal_fw = self.create_autoparal_fw(self.ion_task, ion_fw.spec['_launch_dir'])
        #     self.wf = Workflow([autoparal_fw, ion_fw, ioncell_fw],
        #                        {autoparal_fw: [ion_fw], ion_fw: [ioncell_fw]})
        # else:
        self.wf = Workflow([ion_fw, ioncell_fw], {ion_fw: [ioncell_fw]})

# class RelaxFWWorkflow(AbstractFWWorkflow):
#     """
#     Workflow for structural relaxations. The first task relaxes the atomic position
#     while keeping the unit cell parameters fixed. The second task uses the final
#     structure to perform a structural relaxation in which both the atomic positions
#     and the lattice parameters are optimized.
#     """
#
#     def __init__(self, structure, pseudos, kppa=1000, ecut=None, pawecutdg=None, accuracy="normal",
#                  spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
#                  autoparal=False, folder=None, aim_dilatmx=1.01, **extra_abivars):
#         """
#         Args:
#
#         """
#
#         abiinput = ion_ioncell_relax_input(structure=structure, pseudos=pseudos, kppa=kppa, ecut=ecut,
#                                            pawecutdg=pawecutdg, accuracy=accuracy, spin_mode=spin_mode,
#                                            smearing=smearing, charge=charge, scf_algorithm=scf_algorithm)
#         abiinput.set_vars(**extra_abivars)
#         ion_input, ioncell_input = abiinput.split_datasets()
#
#         abi_ion_task = RelaxTask(StrategyWithInput(ion_input))
#
#         # Create the ionic relaxation fw
#         self.ion_task = RelaxFWTask(abi_ion_task)
#
#         spec = {}
#         if folder:
#             spec['_launch_dir'] = os.path.join(folder, 'atomic_relax')
#
#         ion_fw = Firework(self.ion_task, spec=spec)
#
#         abi_ioncell_task = RelaxTask(StrategyWithInput(ioncell_input))
#
#         # Create the ionic+cell relaxation fw
#         self.ioncell_task = RelaxFWTask(abi_ioncell_task, deps={id(self.ion_task): "structure"})
#
#         spec = {'create_file': 'final'}
#         if folder:
#             spec['_launch_dir'] = os.path.join(folder, 'atomic_and_cell_relax')
#
#         ioncell_fw = Firework(self.ioncell_task, spec=spec)
#
#         # Create the workflow
#         if autoparal:
#             autoparal_fw = self.create_autoparal_fw(self.ion_task, ion_fw.spec['_launch_dir'])
#             self.wf = Workflow([autoparal_fw, ion_fw, ioncell_fw],
#                                {autoparal_fw: [ion_fw], ion_fw: [ioncell_fw]})
#         else:
#             self.wf = Workflow([ion_fw, ioncell_fw], {ion_fw: [ioncell_fw]})


class MultiStepRelaxFWWorkflow(AbstractFWWorkflow):
    """
    Workflow for structural relaxations performed in a multi-step procedure.
    Using a small number of maximum relaxation steps the overall process is
    automatically split in as many FWs are needed.
    """
    def __init__(self, structure, pseudos, ksampling=1000, relax_algo="atoms_only", accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 autoparal=False, max_restart=10, folder=None, **extra_abivars):
        """
        Args:

        """

        task = MultiStepRelaxStrategyFireTask(structure=structure, pseudos=pseudos, ksampling=ksampling,
                                              relax_algo=relax_algo, accuracy=accuracy, spin_mode=spin_mode,
                                              smearing=smearing, charge=charge, scf_algorithm=scf_algorithm,
                                              deps={}, additional_steps=max_restart, **extra_abivars)

        spec = {}
        if folder:
            spec['_launch_dir'] = os.path.join(folder, 'relax')

        fw = Firework(task, spec=spec)

        # Create the workflow
        if autoparal:
            autoparal_fw = self.create_autoparal_fw(task, spec['_launch_dir'])
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
                 folder=None, **extra_abivars):

        # workflow to obtain the density
        if relax_algo in relaxation_methods.keys():
            if relax_algo == "atoms_and_cell":
                dens_wf = RelaxFWWorkflow(structure=structure, pseudos=pseudos, ksampling=ksampling, accuracy=accuracy,
                                          spin_mode=spin_mode, smearing=smearing, charge=charge,
                                          scf_algorithm=scf_algorithm, autoparal=autoparal, folder=folder,
                                          **extra_abivars)
                last_task = dens_wf.ioncell_task
            else:
                dens_wf = RelaxAtomsFWWorkflow(structure=structure, pseudos=pseudos, ksampling=ksampling,
                                               accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
                                               scf_algorithm=scf_algorithm, autoparal=autoparal, folder=folder,
                                               **extra_abivars)
                last_task = dens_wf.task
            deps = {id(last_task): 'DEN structure'}
        else:
            if relax_algo not in (None, 'scf', 'scf_only', 'no_relax'):
                logger.warning('relax_algo value "%s" not recognized. Perorming a scf calculation.' % relax_algo)

            dens_wf = ScfFWWorkflow(structure=structure, pseudos=pseudos, ksampling=ksampling, accuracy=accuracy,
                                    spin_mode=spin_mode, smearing=smearing, charge=charge, scf_algorithm=scf_algorithm,
                                    use_symmetries=use_symmetries, autoparal=autoparal, folder=folder, **extra_abivars)
            last_task = dens_wf.task
            deps = {id(last_task): 'DEN'}

        # Create the ksampling based on high symmetry lines of the structure
        if isinstance(ksampling_band, int):
            # ksampling_band = KSampling.automatic_density(structure, 100)
            ksampling_band = KSampling.path_from_structure(ksampling_band, structure)

        # Create the nscf FW
        nscf_task = NscfStrategyFireTask(scf_task=last_task, ksampling=ksampling_band, nscf_bands=nscf_bands,
                                         nscf_algorithm=nscf_algorithm, deps=deps, **extra_abivars)

        spec = {'create_file': 'band_structure'}
        if folder:
            spec['_launch_dir'] = os.path.join(folder, 'band_structure')

        nscf_fw = Firework(nscf_task, spec=spec)

        # Expand the first part to get the proper list of FWs and links
        fws, links = parse_workflow([dens_wf.wf, nscf_fw], {dens_wf.wf: [nscf_fw]})

        self.wf = Workflow(fws, links)

