# coding: utf-8
"""Factory functions for Abinit input files """
from __future__ import print_function, division, unicode_literals

import numpy as np
import pymatgen.io.abinitio.abiobjects as aobj

from collections import namedtuple
from monty.collections import AttrDict
from pymatgen.io.abinitio.pseudos import PseudoTable
from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput, MultiDataset

import logging
logger = logging.getLogger(__file__)


__all__ = [
    "gs_input",
    "ebands_input",
    "g0w0_with_ppmodel_inputs",
    "bse_with_mdf_inputs",
    "ion_ioncell_relax_input",
    "scf_phonons_inputs",
]


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


# Default values used if user do not specify them
# TODO: Design an object similar to DictVaspInputSet
_DEFAULTS = dict(
    kppa=1000,
)


def _stopping_criterion(runlevel, accuracy):
    """Return the stopping criterion for this runlevel with the given accuracy."""
    tolname = _runl2tolname[runlevel]
    return {tolname: getattr(_tolerances[tolname], accuracy)}


def _find_ecut_pawecutdg(ecut, pawecutdg, pseudos, accuracy='normal'):
    """Return a :class:`AttrDict` with the value of ecut and pawecutdg"""
    # Get ecut and pawecutdg from the pseudo hints.
    has_hints = all(p.has_hints for p in pseudos)

    if ecut is None:
        if has_hints:
            ecut = max(p.hint_for_accuracy(accuracy).ecut for p in pseudos)
        else:
            raise AbiniInput.Error("ecut is None but pseudos do not provide hints for ecut")

    # TODO: This should be the new API.
    if pawecutdg is None and any(p.ispaw for p in pseudos):
        if has_hints:
            pawecutdg = max(p.hint_for_accuracy(accuracy).pawecutdg for p in pseudos)
        else:
            raise RuntimeError("pawecutdg is None but pseudos do not provide hints")

    return AttrDict(ecut=ecut, pawecutdg=pawecutdg)


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
        nband += 10
    else:
        nband += 4

    nband += nband % 2
    return int(nband)


def gs_input(structure, pseudos, 
             kppa=None, ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
             smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
    """
    Returns a :class:`AbinitInput` for band structure calculations.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable` object.
        kppa: Defines the sampling used for the SCF run. Defaults to 1000 if not given.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos
            according to accuracy)
        scf_nband: Number of bands for SCF run. If scf_nband is None, nband is automatically initialized
            from the list of pseudos, the structure and the smearing option.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
    """
    multi = ebands_input(structure, pseudos, 
                 kppa=kppa, 
                 ecut=ecut, pawecutdg=pawecutdg, scf_nband=scf_nband, accuracy=accuracy, spin_mode=spin_mode,
                 smearing=smearing, charge=charge, scf_algorithm=scf_algorithm)

    return multi[0]


def ebands_input(structure, pseudos, 
                 kppa=None, nscf_nband=None, ndivsm=15, 
                 ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
                 smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, dos_kppa=None):
    """
    Returns a :class:`AbinitInput` for band structure calculations.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable` object.
        kppa: Defines the sampling used for the SCF run. Defaults to 1000 if not given.
        nscf_nband: Number of bands included in the NSCF run. Set to scf_nband + 10 if None.
        ndivsm: Number of divisions used to sample the smallest segment of the k-path.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos
            according to accuracy)
        scf_nband: Number of bands for SCF run. If scf_nband is None, nband is automatically initialized
            from the list of pseudos, the structure and the smearing option.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
        dos_kppa: Scalar or List of integers with the number of k-points per atom
            to be used for the computation of the DOS (None if DOS is not wanted).
    """
    structure = Structure.as_structure(structure)

    if dos_kppa is not None and not isinstance(dos_kppa, (list, tuple)):
        dos_kppa = [dos_kppa]

    multi = MultiDataset(structure, pseudos, ndtset=2 if dos_kppa is None else 2 + len(dos_kppa))

    # Set the cutoff energies.
    multi.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos))

    # SCF calculation.
    kppa = _DEFAULTS.get("kppa") if kppa is None else kppa
    scf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm, 
                                   charge=charge, nband=scf_nband, fband=None)

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, multi.pseudos, scf_electrons)

    multi[0].set_vars(scf_ksampling.to_abivars())
    multi[0].set_vars(scf_electrons.to_abivars())
    multi[0].set_vars(_stopping_criterion("scf", accuracy))

    # Band structure calculation.
    nscf_ksampling = aobj.KSampling.path_from_structure(ndivsm, structure)
    nscf_nband = scf_electrons.nband + 10 if nscf_nband is None else nscf_nband
    nscf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
                                    charge=charge, nband=nscf_nband, fband=None)

    multi[1].set_vars(nscf_ksampling.to_abivars())
    multi[1].set_vars(nscf_electrons.to_abivars())
    multi[1].set_vars(_stopping_criterion("nscf", accuracy))

    # DOS calculation with different values of kppa.
    if dos_kppa is not None:
        for i, kppa in enumerate(dos_kppa):
            dos_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
            #dos_ksampling = aobj.KSampling.monkhorst(dos_ngkpt, shiftk=dos_shiftk, chksymbreak=0)
            dos_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
                                           charge=charge, nband=nscf_nband) 
            dt = 2 + i
            multi[dt].set_vars(dos_ksampling.to_abivars())
            multi[dt].set_vars(dos_electrons.to_abivars())
            multi[dt].set_vars(_stopping_criterion("nscf", accuracy))

    return multi


def ion_ioncell_relax_input(structure, pseudos, 
                            kppa=None, nband=None,
                            ecut=None, pawecutdg=None, accuracy="normal", spin_mode="polarized",
                            smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
    """
    Returns a :class:`AbinitInput` for a structural relaxation. The first dataset optmizes the 
    atomic positions at fixed unit cell. The second datasets optimizes both ions and unit cell parameters.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable` object.
        kppa: Defines the sampling used for the Brillouin zone.
        nband: Number of bands included in the SCF run.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
    """
    structure = Structure.as_structure(structure)
    multi = MultiDataset(structure, pseudos, ndtset=2)

    # Set the cutoff energies.
    multi.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos))

    kppa = _DEFAULTS.get("kppa") if kppa is None else kppa
    ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
    electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm, 
                               charge=charge, nband=nband, fband=None)

    if electrons.nband is None:
        electrons.nband = _find_scf_nband(structure, multi.pseudos, electrons)

    ion_relax = aobj.RelaxationMethod.atoms_only(atoms_constraints=None)
    ioncell_relax = aobj.RelaxationMethod.atoms_and_cell(atoms_constraints=None)

    multi.set_vars(electrons.to_abivars())
    multi.set_vars(ksampling.to_abivars())

    multi[0].set_vars(ion_relax.to_abivars())
    multi[0].set_vars(_stopping_criterion("relax", accuracy))

    multi[1].set_vars(ioncell_relax.to_abivars())
    multi[1].set_vars(_stopping_criterion("relax", accuracy))

    return multi


def ion_ioncell_relax_and_ebands_input(structure, pseudos, 
                                        kppa=None, nband=None,
                                        ecut=None, pawecutdg=None, accuracy="normal", spin_mode="polarized",
                                        smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
    """
    Returns a :class:`AbinitInput` for a structural relaxation. The first dataset optmizes the 
    atomic positions at fixed unit cell. The second datasets optimizes both ions and unit cell parameters.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable` object.
        kppa: Defines the sampling used for the Brillouin zone.
        nband: Number of bands included in the SCF run.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
    """
    structure = Structure.as_structure(structure)

    relax_multi = ion_ioncell_relax_input(structure, pseudos, 
                                          kppa=kppa, nband=nband,
                                          ecut=ecut, pawecutdg=pawecutdg, accuracy=accuracy, spin_mode=spin_mode,
                                          smearing=smearing, charge=charge, scf_algorithm=scf_algorithm)

    ebands_multi = ebands_input(structure, pseudos, 
                                kppa=kppa, nscf_nband=None, ndivsm=15, 
                                ecut=ecut, pawecutdg=pawecutdg, scf_nband=None, accuracy=accuracy, spin_mode=spin_mode,
                                smearing=smearing, charge=charge, scf_algorithm=scf_algorithm, dos_kppa=None)

    return relax_multi + ebands_multi


def g0w0_with_ppmodel_inputs(structure, pseudos, 
                            kppa, nscf_nband, ecuteps, ecutsigx,
                            ecut=None, pawecutdg=None,
                            accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
                            ppmodel="godby", charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None,
                            sigma_nband=None, gw_qprange=1):
    """
    Returns a :class:`AbinitInput` object that performs G0W0 calculations with the plasmon pole approximation.

    Args:
        structure: Pymatgen structure.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable` object.
        kppa: Defines the sampling used for the SCF run.
        nscf_nband: Number of bands included in the NSCF run.
        ecuteps: Cutoff energy [Ha] for the screening matrix.
        ecutsigx: Cutoff energy [Ha] for the exchange part of the self-energy.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized 
            from the pseudos according to accuracy)
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
    multi = MultiDataset(structure, pseudos, ndtset=4)

    # Set the cutoff energies.
    multi.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos))

    scf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm, 
                                   charge=charge, nband=None, fband=None)

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, multi.pseudos, scf_electrons)

    multi[0].set_vars(scf_ksampling.to_abivars())
    multi[0].set_vars(scf_electrons.to_abivars())
    multi[0].set_vars(_stopping_criterion("scf", accuracy))

    nscf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
    nscf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
                                    charge=charge, nband=nscf_nband, fband=None)

    multi[1].set_vars(nscf_ksampling.to_abivars())
    multi[1].set_vars(nscf_electrons.to_abivars())
    multi[1].set_vars(_stopping_criterion("nscf", accuracy))
    # nbdbuf

    # Screening.
    if scr_nband is None: scr_nband = nscf_nband
    screening = aobj.Screening(ecuteps, scr_nband, w_type="RPA", sc_mode="one_shot",
                               hilbert=None, ecutwfn=None, inclvkb=inclvkb)

    multi[2].set_vars(nscf_ksampling.to_abivars())
    multi[2].set_vars(nscf_electrons.to_abivars())
    multi[2].set_vars(screening.to_abivars())
    multi[2].set_vars(_stopping_criterion("screening", accuracy)) # Dummy
    #scr_strategy = ScreeningStrategy(scf_strategy, nscf_strategy, screening)

    # Sigma.
    if sigma_nband is None: sigma_nband = nscf_nband
    self_energy = aobj.SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening,
                             gw_qprange=gw_qprange, ppmodel=ppmodel)

    multi[3].set_vars(nscf_ksampling.to_abivars())
    multi[3].set_vars(nscf_electrons.to_abivars())
    multi[3].set_vars(self_energy.to_abivars())
    multi[3].set_vars(_stopping_criterion("sigma", accuracy)) # Dummy
    #sigma_strategy = aobj.SelfEnergyStrategy(scf_strategy, nscf_strategy, scr_strategy, self_energy)

    # TODO: Cannot use istwfk != 1.
    multi.set_vars(istwfk="*1")

    return multi

#TODO
#def g0w0_extended_work(structure, pseudos, kppa, nscf_nband, ecuteps, ecutsigx, scf_nband, accuracy="normal",

def bse_with_mdf_inputs(structure, pseudos, 
                        scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk, 
                        ecuteps, bs_loband, bs_nband, soenergy, mdf_epsinf, 
                        ecut=None, pawecutdg=None, 
                        exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="polarized", 
                        smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
    """
    Returns a :class:`AbinitInput` object that performs a GS + NSCF + Bethe-Salpeter calculation.
    The self-energy corrections are approximated with the scissors operator.
    The screening in modeled with the model dielectric function.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable` object.
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
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos
            according to accuracy)
        exc_type: Approximation used for the BSE Hamiltonian (Tamm-Dancoff or coupling).
        bs_algo: Algorith for the computatio of the macroscopic dielectric function.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving the SCF cycle.
    """
    structure = Structure.as_structure(structure)
    multi = MultiDataset(structure, pseudos, ndtset=3)

    # Set the cutoff energies.
    d = _find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos)
    multi.set_vars(ecut=d.ecut, ecutwfn=d.ecut, pawecutdg=d.pawecutdg)

    # Ground-state 
    scf_ksampling = aobj.KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)

    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm, 
                                   charge=charge, nband=None, fband=None)

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, multi.pseudos, scf_electrons)

    multi[0].set_vars(scf_ksampling.to_abivars())
    multi[0].set_vars(scf_electrons.to_abivars())
    multi[0].set_vars(_stopping_criterion("scf", accuracy))

    # NSCF calculation with the randomly-shifted k-mesh.
    nscf_ksampling = aobj.KSampling.monkhorst(nscf_ngkpt, shiftk=nscf_shiftk, chksymbreak=0)

    nscf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
                                    charge=charge, nband=nscf_nband, fband=None)

    multi[1].set_vars(nscf_ksampling.to_abivars())
    multi[1].set_vars(nscf_electrons.to_abivars())
    multi[1].set_vars(_stopping_criterion("nscf", accuracy))

    # BSE calculation.
    exc_ham = aobj.ExcHamiltonian(bs_loband, bs_nband, soenergy, coulomb_mode="model_df", ecuteps=ecuteps, 
                                  spin_mode=spin_mode, mdf_epsinf=mdf_epsinf, exc_type=exc_type, algo=bs_algo,
                                  bs_freq_mesh=None, with_lf=True, zcut=None)

    multi[2].set_vars(nscf_ksampling.to_abivars())
    multi[2].set_vars(nscf_electrons.to_abivars())
    multi[2].set_vars(exc_ham.to_abivars())
    #multi[2].set_vars(_stopping_criterion("nscf", accuracy))

    # TODO: Cannot use istwfk != 1.
    multi.set_vars(istwfk="*1")

    return multi


def scf_phonons_inputs(structure, pseudos, kppa,
                       ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
                       smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):

    """
    Returns a :class:`AbinitInput` for performing phonon calculations.
    GS input + the input files for the phonon calculation.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable` object.
        kppa: Defines the sampling used for the SCF run.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the
            pseudos according to accuracy)
        scf_nband: Number of bands for SCF run. If scf_nband is None, nband is automatically initialized from the list of 
            pseudos, the structure and the smearing option.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
    """
    # Build the input file for the GS run.
    gs_inp = AbinitInput(structure=structure, pseudos=pseudos)

    # Set the cutoff energies.
    gs_inp.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, gs_inp.pseudos))

    ksampling = aobj.KSampling.automatic_density(gs_inp.structure, kppa, chksymbreak=0)
    gs_inp.set_vars(ksampling.to_abivars())
    gs_inp.set_vars(tolvrs=1.0e-18)

    # Get the qpoints in the IBZ. Note that here we use a q-mesh with ngkpt=(4,4,4) and shiftk=(0,0,0)
    # i.e. the same parameters used for the k-mesh in gs_inp.
    qpoints = gs_inp.abiget_ibz(ngkpt=(4,4,4), shiftk=(0,0,0), kptopt=1).points
    #print("get_ibz qpoints:", qpoints)

    # Build the input files for the q-points in the IBZ.
    ph_inputs = MultiDataset(gs_inp.structure, pseudos=gs_inp.pseudos, ndtset=len(qpoints))

    for ph_inp, qpt in zip(ph_inputs, qpoints):
        # Response-function calculation for phonons.
        ph_inp.set_vars(
            rfphon=1,        # Will consider phonon-type perturbation
            nqpt=1,          # One wavevector is to be considered
            qpt=qpt,         # This wavevector is q=0 (Gamma)
            tolwfr=1.0e-20,
            kptopt=3,        # One could used symmetries for Gamma.
        )
            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            #kptopt   2      # Automatic generation of k points, taking

        irred_perts = ph_inp.abiget_irred_phperts()

        #for pert in irred_perts:
        #    #print(pert)
        #    # TODO this will work for phonons, but not for the other types of perturbations.
        #    ph_inp = q_inp.deepcopy()

        #    rfdir = 3 * [0]
        #    rfdir[pert.idir -1] = 1

        #    ph_inp.set_vars(
        #        rfdir=rfdir,
        #        rfatpol=[pert.ipert, pert.ipert]
        #    )

        #    ph_inputs.append(ph_inp)

    # Split input into gs_inp and ph_inputs
    all_inps = [gs_inp] 
    all_inps.extend(ph_inputs.split_datasets())

    return all_inps


def phonons_from_gsinput(gs_inp):
    """
    Returns a :class:`AbinitInput` for performing phonon calculations.
    GS input + the input files for the phonon calculation.
    """
    qpoints = gs_inp.abiget_ibz(ngkpt=(4,4,4), shiftk=(0,0,0), kptopt=1).points
    #print("get_ibz qpoints:", qpoints)

    # Build the input files for the q-points in the IBZ.
    # Response-function calculation for phonons.
    ph_inputs = []
    for qpt in qpoints:
        q_inp = gs_inp.deepcopy()
        #q_inp.pop_tolerances()
        q_inp.set_vars(
            rfphon=1,        # Will consider phonon-type perturbation
            nqpt=1,          # One wavevector is to be considered
            qpt=qpt,         # This wavevector is q=0 (Gamma)
            #tolwfr=1.0e-20,
            kptopt=3,        # One could used symmetries for Gamma.
        )
            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            #kptopt   2      # Automatic generation of k points, taking
        #print(tmp_inp)
        irred_perts = q_inp.abiget_irred_phperts()

        for pert in irred_perts:
            #print(pert)
            # TODO this will work for phonons, but not for the other types of perturbations.
            ph_inp = q_inp.deepcopy()

            rfdir = 3 * [0]
            rfdir[pert.idir -1] = 1

            ph_inp.set_vars(
                rfdir=rfdir,
                rfatpol=[pert.ipert, pert.ipert]
            )

            ph_inputs.append(ph_inp)

    #return MultiDataset.from_inputs(ph_inputs)
