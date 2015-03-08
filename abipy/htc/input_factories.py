# coding: utf-8
"""Factory function for Abinit input files """
from __future__ import print_function, division, unicode_literals

import six
import abc
import numpy as np
import pymatgen.io.abinitio.abiobjects as aobj 

from collections import OrderedDict, namedtuple 
from monty.string import is_string, list_strings
from monty.collections import dict2namedtuple
from pymatgen.io.abinitio.pseudos import PseudoTable
#from pymatgen.serializers.json_coders import PMGSONable
from abipy.core.structure import Structure
from .input import AbiInput

import logging
logger = logging.getLogger(__file__)


# TODO: To be discussed: 
#    1) extra_abivars is more similar to a hack. The factory functions are designed for
#       HPC hence we cannot allow the user to inject something we cannot control easily
#       Shall we remove it?
#    2) scf_nband and nscf_band should be computed from the pseudos, the structure
#       and some approximation for the band dispersion.
#       SCF fails if nband is too small or has problems if we don't have enough partially
#       occupied states in metals (can write EventHandler but it would be nice if we could
#       fix this problem in advance.
#    3) How do we handle options related to parallelism e.g. paral_kgb?
#    4) The API of the factory functions must be simple enough so that we can easily generate
#       flows but, on the other hand, we would like to decorate the input with extra features
#       e.g. we would like to do a LDA+U band structure, a LDA+U relaxation etc.
#       For a possible solution based on factory functions see:
#
#            http://python-3-patterns-idioms-test.readthedocs.org/en/latest/Factory.html
#
#       for decorator pattern see:
#
#            http://www.tutorialspoint.com/design_pattern/decorator_pattern.htm

class DecoratorError(Exception):
    """Error class raised by :class:`InputDecorator`."""


class InputDecorator(six.with_metaclass(abc.ABCMeta, object)):
    """Abstract Base class."""
    Error = DecoratorError

    @abc.abstractmethod
    def as_dict(self):

    def decorate(self, inp, deepcopy=True):
        new_inp = self._decorate(inp, deepcopy=deepcopy)
        # Log the decoration in new_inp.
        new_inp._decorators.append(self.as_dict())
        return new_inp

    @abc.abstractmethod
    def _decorate(self, inp, deepcopy=True):
        """

        Args:
            inp: :class:`AbiInput` object.
            deepcopy: True if a deepcopy of inp should be performed 
                before changing the object.

        Returns:
            decorated :class:`AbiInput` object (new object)
        """

#class SpinDecorator(InputDecorator):
#    def __init__(self, spinmode):
#        """Change the spin polarization."""
#        self.spinmode = aobj.SpinMpde.as_spinmode(spin_mode)
#
#    def _decorate(self, inp, deepcopy=True)
#        if deepcopy: inp = inp.deepcopy()
#        inp.set_vars(self.spinmode.to_abivars())
#        return inp

#class SmearingDecorator(InputDecorator):
#    """Change the electronic smearing."""
#    def __init__(self, spinmode):
#        self.smearing = aobj.Smearing.as_smearing(smearing)
#
#    def _decorate(self, inp, deepcopy=True)
#        if deepcopy: inp = inp.deepcopy()
#        inp.set_vars(self.smearing.to_abivars())
#        return inp

#class UJDecorator(InputDecorator):
#    """Add LDA+U to an :class:`AbiInput` object."""
#    def __init__(self, luj_for_symbol, usepawu=1):
#        self.usepawu = usepawu
#        self.luj_for_symbol = luj_for_symbol
#
#    def _decorate(self, inp, deepcopy=True)
#        if not inp.ispaw: raise self.Error("LDA+U requires PAW!")
#        if deepcopy: inp = inp.deepcopy()
#
#        inp.set_vars(usepawu=self.usepawu)
#        return inp


#class LexxDecorator(InputDecorator):
#    """Add LDA+U to an :class:`AbiInput` object."""
#    def __init__(self, luj_for_symbol, usepawu=1):
#        self.usepawu = usepawu
#        self.luj_for_symbol = luj_for_symbol
#
#    def _decorate(self, inp, deepcopy=True)
#        if not inp.ispaw: raise self.Error("LDA+U requires PAW!")
#        if deepcopy: inp = inp.deepcopy()
#
#        inp.set_vars(usepawu=self.usepawu)
#        return inp


# Name of the (default) tolerance used by the runlevels.
_runl2tolname = {
    "scf": 'tolvrs',
    "nscf": 'tolwfr',
    "dfpt": 'toldfe',        # ?
    "screening": 'toldfe',   # dummy
    "sigma": 'toldfe',       # dummy
    "bse": 'toldfe',         # ?
    "relax": 'tolrff',
}

# Tolerances for the different levels of accuracy.
T = namedtuple('Tolerance', "low normal high")
_tolerances = {
    "toldfe": T(1.e-7,  1.e-8,  1.e-9),
    "tolvrs": T(1.e-7,  1.e-8,  1.e-9),
    "tolwfr": T(1.e-15, 1.e-17, 1.e-19),
    "tolrff": T(0.04,   0.02,   0.01)}
del T


def _stopping_criterion(runlevel, accuracy):
    """Return the stopping criterion for this runlevel with the given accuracy."""
    tolname = _runl2tolname[runlevel]
    return {tolname: getattr(_tolerances[tolname], accuracy)}


def _find_ecut_pawecutdg(ecut, pawecutdg, pseudos):
    """Return the value of ecut and pawecutdg"""
    #TODO: Get ecut and pawecutdg from the pseudo hints.
    #if ecut is None: 
    #  ecut = max(p.hint for p in pseudos)

    #if pawecutdg is None and any(p.ispaw for p in pseudos):
    #    pawecutdg = max(p.hint for p in pseudos)

    return ecut, pawecutdg


def _find_scf_nband(structure, pseudos, electrons):
    """Find the value of nband."""
    if electrons.nband is not None: return electrons.nband

    nsppol, smearing = electrons.nsppol, electrons.smearing
    
    # Number of valence electrons including possible extra charge
    nval = structure.num_valence_electrons(pseudos)
    nval -= electrons.charge

    # First guess (semiconductors)
    nband = nval // nsppol

    # TODO: Find better algorithm
    # If nband is too small we may kill the job, increase nband and restart
    # but this change could cause problems in the other steps of the calculation
    # if the change is not propagated e.g. phonons in metals.
    if smearing:
        # metallic occupation
        nband += 12
    else:
        nband += 4

    nband += nband % 2
    return int(nband)


def ebands_input(structure, pseudos, scf_kppa, nscf_nband, ndivsm, 
                 ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
                 smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, dos_kppa=None):
    """
    Returns a :class:`AbiInput` for band structure calculations.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable: object.
        scf_kppa: Defines the sampling used for the SCF run.
        nscf_nband: Number of bands included in the NSCF run.
        ndivsm: Number of divisions used to sample the smallest segment of the k-path.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos according to accuracy)
        scf_nband: Number of bands for SCF run. If scf_nband is None, nband is automatically initialized from the list of 
            pseudos, the structure and the smearing option.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
        dos_kppa: Scalar or List of integers with the number of k-points per atom
            to be used for the computation of the DOS (None if DOS is not wanted).
    """
    structure = Structure.as_structure(structure)
    pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(structure)

    if dos_kppa is not None and not isinstance(dos_kppa, (list, tuple)):
        dos_kppa = [dos_kppa]

    inp = AbiInput(pseudos, ndtset=2 if dos_kppa is None else 2 + len(dos_kppa))
    inp.set_structure(structure)

    # Set the cutoff energy.
    ecut, pawecutdg = _find_ecut_pawecutdg(ecut, pawecutdg, pseudos)
    inp.set_vars(ecut=ecut, pawecutdg=pawecutdg)

    # SCF calculation.
    scf_ksampling = aobj.KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm, 
                                   charge=charge, nband=scf_nband, fband=None)

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, pseudos, scf_electrons)

    inp[1].set_vars(scf_ksampling.to_abivars())
    inp[1].set_vars(scf_electrons.to_abivars())
    inp[1].set_vars(_stopping_criterion("scf", accuracy))

    # Band structure calculation.
    nscf_ksampling = aobj.KSampling.path_from_structure(ndivsm, structure)
    nscf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
                                    charge=charge, nband=nscf_nband, fband=None)

    inp[2].set_vars(nscf_ksampling.to_abivars())
    inp[2].set_vars(nscf_electrons.to_abivars())
    inp[2].set_vars(_stopping_criterion("nscf", accuracy))

    # DOS calculation with different values of kppa.
    if dos_kppa is not None:
        for i, kppa in enumerate(dos_kppa):
            dos_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
            #dos_ksampling = aobj.KSampling.monkhorst(dos_ngkpt, shiftk=dos_shiftk, chksymbreak=0)
            dos_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
                                           charge=charge, nband=nscf_nband) 
            dt = 3 + i
            inp[dt].set_vars(dos_ksampling.to_abivars())
            inp[dt].set_vars(dos_electrons.to_abivars())
            inp[dt].set_vars(_stopping_criterion("nscf", accuracy))

    return inp
    #scf_inp, nscf_inp, dos_inps = inp.split_datasets()
    #return dict2namedtuple(scf_inp=scf_inp, nscf_inp=nscf_inp, dos_inps")


def ion_ioncell_relax_input(structure, pseudos, kppa, nband=None,
                            ecut=None, pawecutdg=None, accuracy="normal", spin_mode="polarized",
                            smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
    """
    Returns a :class:`AbiInput` for a structural relaxation. The first dataset optmizes the 
    atomic positions at fixed unit cell. The second datasets optimizes both ions and unit cell parameters.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable: object.
        kppa: Defines the sampling used for the Brillouin zone.
        nband: Number of bands included in the SCF run.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
    """
    structure = Structure.as_structure(structure)
    pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(structure)

    inp = AbiInput(pseudos, ndtset=2)
    inp.set_structure(structure)

    # Set the cutoff energy
    ecut, pawecutdg = _find_ecut_pawecutdg(ecut, pawecutdg, pseudos)
    inp.set_vars(ecut=ecut, pawecutdg=pawecutdg)

    ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
    electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm, 
                               charge=charge, nband=nband, fband=None)

    if electrons.nband is None:
        electrons.nband = _find_scf_nband(structure, pseudos, electrons)

    ion_relax = aobj.RelaxationMethod.atoms_only(atoms_constraints=None)
    ioncell_relax = aobj.RelaxationMethod.atoms_and_cell(atoms_constraints=None)

    inp.set_vars(electrons.to_abivars())
    inp.set_vars(ksampling.to_abivars())

    inp[1].set_vars(ion_relax.to_abivars())
    # Stopping criterion is already in RelaxationMethod
    # TODO: Use similar approach for Scf
    #inp[1].set_vars(_stopping_criterion("relax", accuracy))

    inp[2].set_vars(ioncell_relax.to_abivars())

    return inp
    #return dict2namedtuple(scf_inp=scf_inp, nscf_inp=nscf_inp, dos_inps")


def g0w0_with_ppmodel_input(structure, pseudos, scf_kppa, nscf_nband, ecuteps, ecutsigx,
                            ecut=None, pawecutdg=None,
                            accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
                            ppmodel="godby", charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None,
                            sigma_nband=None, gw_qprange=1):
    """
    Returns a :class:`AbiInput` object that performs G0W0 calculations with the plasmon pole approximation.

    Args:
        structure: Pymatgen structure.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable: object.
        scf_kppa: Defines the sampling used for the SCF run.
        nscf_nband: Number of bands included in the NSCF run.
        ecuteps: Cutoff energy [Ha] for the screening matrix.
        ecutsigx: Cutoff energy [Ha] for the exchange part of the self-energy.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos according to accuracy)
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        ppmodel: Plasmonpole technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
        inclvkb: Treatment of the dipole matrix elements (see abinit variable).
        scr_nband: Number of bands used to compute the screening (default is nscf_nband)
        sigma_nband: Number of bands used to compute the self-energy (default is nscf_nband)
        gw_qprange: Option for the automatic selection of k-points and bands for GW corrections.
            See Abinit docs for more detail. The default value makes the code compute the
            QP energies for all the point in the IBZ and one band above and one band below the Fermi level.
    """
    structure = Structure.as_structure(structure)
    pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(structure)

    inp = AbiInput(pseudos, ndtset=4)
    inp.set_structure(structure)

    # Set the cutoff energy
    ecut, pawecutdg = _find_ecut_pawecutdg(ecut, pawecutdg, pseudos)
    inp.set_vars(ecut=ecut, pawecutdg=pawecutdg)

    scf_ksampling = aobj.KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm, 
                                   charge=charge, nband=None, fband=None)

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, pseudos, scf_electrons)

    inp[1].set_vars(scf_ksampling.to_abivars())
    inp[1].set_vars(scf_electrons.to_abivars())
    inp[1].set_vars(_stopping_criterion("scf", accuracy))

    nscf_ksampling = aobj.KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
    nscf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
                                    charge=charge, nband=nscf_nband, fband=None)

    inp[2].set_vars(nscf_ksampling.to_abivars())
    inp[2].set_vars(nscf_electrons.to_abivars())
    inp[2].set_vars(_stopping_criterion("nscf", accuracy))
    # nbdbuf

    # Screening.
    if scr_nband is None: scr_nband = nscf_nband
    screening = aobj.Screening(ecuteps, scr_nband, w_type="RPA", sc_mode="one_shot",
                          hilbert=None, ecutwfn=None, inclvkb=inclvkb)
    inp[3].set_vars(screening.to_abivars())
    #scr_strategy = ScreeningStrategy(scf_strategy, nscf_strategy, screening)

    # Sigma.
    if sigma_nband is None: sigma_nband = nscf_nband
    self_energy = aobj.SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening,
                             gw_qprange=gw_qprange, ppmodel=ppmodel)
    inp[4].set_vars(self_energy.to_abivars())
    #sigma_strategy = aobj.SelfEnergyStrategy(scf_strategy, nscf_strategy, scr_strategy, self_energy)

    # TODO: Cannot use istwfk != 1.
    inp.set_vars(istwfk="*1")

    return inp
    #return dict2namedtuple(scf_inp=scf_inp, nscf_inp=nscf_inp, dos_inps")

#TODO
#def g0w0_extended_work(structure, pseudos, scf_kppa, nscf_nband, ecuteps, ecutsigx, scf_nband, accuracy="normal",


def bse_with_mdf_input(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk, 
                       ecuteps, bs_loband, bs_nband, soenergy, mdf_epsinf, 
                       ecut=None, pawecutdg=None, 
                       exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="polarized", 
                       smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
    """
    Returns a :class:`AbiInput` object that performs a GS + NSCF + Bethe-Salpeter calculation.
    The self-energy corrections are approximated with the scissors operator.
    The screening in modeled with the model dielectric function.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable: object.
        scf_kppa: Defines the sampling used for the SCF run.
        nscf_nband: Number of bands included in the NSCF run.
        nscf_ngkpt: Divisions of the k-mesh used for the NSCF and the BSE run.
        nscf_shiftk: Shifts used for the NSCF and the BSE run.
        ecuteps: Cutoff energy [Ha] for the screening matrix.
        bs_loband: Index of the first occupied band included the e-h basis set
            (ABINIT convention i.e. first band starts at 1).
            Can be scalar or array of shape (nsppol,)
        bs_nband: Highest band idex used for the construction of the e-h basis set.
        soenergy: Scissor energy in Hartree.
        mdf_epsinf: Value of the macroscopic dielectric function used in expression for the model dielectric function.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos according to accuracy)
        exc_type: Approximation used for the BSE Hamiltonian (Tamm-Dancoff or coupling).
        bs_algo: Algorith for the computatio of the macroscopic dielectric function.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving the SCF cycle.
    """
    structure = Structure.as_structure(structure)
    pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(structure)

    inp = AbiInput(pseudos, ndtset=3)
    inp.set_structure(structure)

    # Set the cutoff energy
    ecut, pawecutdg = _find_ecut_pawecutdg(ecut, pawecutdg, pseudos)
    inp.set_vars(ecut=ecut, pawecutdg=pawecutdg)

    # Ground-state 
    scf_ksampling = aobj.KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)

    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm, 
                                   charge=charge, nband=None, fband=None)

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, pseudos, scf_electrons)

    inp[1].set_vars(scf_ksampling.to_abivars())
    inp[1].set_vars(scf_electrons.to_abivars())
    inp[1].set_vars(_stopping_criterion("scf", accuracy))

    # NSCF calculation with the randomly-shifted k-mesh.
    nscf_ksampling = aobj.KSampling.monkhorst(nscf_ngkpt, shiftk=nscf_shiftk, chksymbreak=0)

    nscf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
                                    charge=charge, nband=nscf_nband, fband=None)

    inp[2].set_vars(nscf_ksampling.to_abivars())
    inp[2].set_vars(nscf_electrons.to_abivars())
    inp[2].set_vars(_stopping_criterion("nscf", accuracy))

    # BSE calculation.
    exc_ham = aobj.ExcHamiltonian(bs_loband, bs_nband, soenergy, coulomb_mode="model_df", ecuteps=ecuteps, 
                                  spin_mode=spin_mode, mdf_epsinf=mdf_epsinf, exc_type=exc_type, algo=bs_algo,
                                  bs_freq_mesh=None, with_lf=True, zcut=None)

    inp[3].set_vars(nscf_ksampling.to_abivars())
    inp[3].set_vars(nscf_electrons.to_abivars())
    inp[3].set_vars(exc_ham.to_abivars())
    #inp[3].set_vars(_stopping_criterion("nscf", accuracy))

    # TODO: Cannot use istwfk != 1.
    inp.set_vars(istwfk="*1")

    return inp
    #return dict2namedtuple(scf_inp=scf_inp, nscf_inp=nscf_inp, dos_inps")


def scf_phonons_inputs(structure, pseudos, scf_kppa,
                       ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
                       smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):

    """
    Returns a :class:`AbiInput` for performing phonon calculations.
    GS input + the input files for the phonon calculation.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable: object.
        scf_kppa: Defines the sampling used for the SCF run.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos according to accuracy)
        scf_nband: Number of bands for SCF run. If scf_nband is None, nband is automatically initialized from the list of 
            pseudos, the structure and the smearing option.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
    """
    structure = Structure.as_structure(structure)
    pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(structure)

    # List of q-points for the phonon calculation.
    #qpoints = np.reshape([
    #         0.00000000E+00,  0.00000000E+00,  0.00000000E+00, 
    #         2.50000000E-01,  0.00000000E+00,  0.00000000E+00,
    #         5.00000000E-01,  0.00000000E+00,  0.00000000E+00,
    #         2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
    #         5.00000000E-01,  2.50000000E-01,  0.00000000E+00,
    #        -2.50000000E-01,  2.50000000E-01,  0.00000000E+00,
    #         5.00000000E-01,  5.00000000E-01,  0.00000000E+00,
    #        -2.50000000E-01,  5.00000000E-01,  2.50000000E-01,
    #        ], (-1, 3))

    # Global variables used both for the GS and the DFPT run.
    #global_vars = dict(nband=4,             
    #                   #ecut=3.0,         
    #                   #ecut=12.0,
    #                   ngkpt=[4, 4, 4],
    #                   nshiftk=4,
    #                   shiftk=[0.0, 0.0, 0.5,   # This gives the usual fcc Monkhorst-Pack grid
    #                           0.0, 0.5, 0.0,
    #                           0.5, 0.0, 0.0,
    #                           0.5, 0.5, 0.5],
    #                   #shiftk=[0, 0, 0],
    #                   #paral_kgb=paral_kgb,
    #                   #ixc=1,
    #                   #nstep=25,
    #                   #diemac=9.0,
    #                )

    # Build the input file for the GS run.
    gs_inp = AbiInput(pseudos=pseudos)
    gs_inp.set_structure(structure)

    # Set the cutoff energy
    ecut, pawecutdg = _find_ecut_pawecutdg(ecut, pawecutdg, pseudos)
    gs_inp.set_vars(ecut=ecut, pawecutdg=pawecutdg)

    ksampling = aobj.KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
    gs_inp.set_vars(ksampling.to_abivars())
    gs_inp.set_vars(tolvrs=1.0e-18)
    #gs_inp.set_vars(global_vars)

    # Get the qpoints in the IBZ. Note that here we use a q-mesh with ngkpt=(4,4,4) and shiftk=(0,0,0)
    # i.e. the same parameters used for the k-mesh in gs_inp.
    qpoints = gs_inp.get_ibz(ngkpt=(4,4,4), shiftk=(0,0,0), kptopt=1).points
    print("get_ibz qpoints:", qpoints)

    # Build the input files for the q-points in the IBZ.
    ph_inputs = AbiInput(pseudos=pseudos, ndtset=len(qpoints))
    ph_inputs.set_structure(structure)

    #pprint(inp[1].get_irred_perts())
    for ph_inp, qpt in zip(ph_inputs, qpoints):
        # Response-function calculation for phonons.
        #ph_inp.set_vars(global_vars)
        ph_inp.set_vars(
            rfphon=1,        # Will consider phonon-type perturbation
            nqpt=1,          # One wavevector is to be considered
            qpt=qpt,         # This wavevector is q=0 (Gamma)
            tolwfr=1.0e-20,
            kptopt=3,
            )
            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            #kptopt   2      # Automatic generation of k points, taking

    # Split input into gs_inp and ph_inputs
    all_inps = [gs_inp] 
    all_inps.extend(ph_inputs.split_datasets())

    return all_inps
    #return dict2namedtuple(scf_inp=scf_inp, nscf_inp=nscf_inp, dos_inps")
