# coding: utf-8
"""Factory functions for Abinit input files """
from __future__ import print_function, division, unicode_literals

import numpy as np
import pymatgen.io.abinit.abiobjects as aobj

from collections import namedtuple
from monty.collections import AttrDict
from monty.json import jsanitize, MontyDecoder
from pymatgen.io.abinit.pseudos import PseudoTable
from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput, MultiDataset
from abipy.abio.input_tags import *

import logging
from pymatgen.serializers.json_coders import pmg_serialize

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
    if ecut is None or (pawecutdg is None and any(p.ispaw for p in pseudos)):
        has_hints = all(p.has_hints for p in pseudos)

    if ecut is None:
        if has_hints:
            ecut = max(p.hint_for_accuracy(accuracy).ecut for p in pseudos)
        else:
            raise AbinitInput.Error("ecut is None but pseudos do not provide hints for ecut")

    # TODO: This should be the new API.
    if pawecutdg is None and any(p.ispaw for p in pseudos):
        if has_hints:
            pawecutdg = max(p.hint_for_accuracy(accuracy).pawecutdg for p in pseudos)
        else:
            raise RuntimeError("pawecutdg is None but pseudos do not provide hints")

    return AttrDict(ecut=ecut, pawecutdg=pawecutdg)


def _find_scf_nband(structure, pseudos, electrons, spinat=None):
    """Find the value of nband."""
    if electrons.nband is not None: return electrons.nband

    nsppol, smearing = electrons.nsppol, electrons.smearing
    
    # Number of valence electrons including possible extra charge
    nval = structure.num_valence_electrons(pseudos)
    nval -= electrons.charge

    # First guess (semiconductors)
    nband = nval // 2

    # TODO: Find better algorithm
    # If nband is too small we may kill the job, increase nband and restart
    # but this change could cause problems in the other steps of the calculation
    # if the change is not propagated e.g. phonons in metals.
    if smearing:
        # metallic occupation
        nband = max(np.ceil(nband*1.2), nband+10)
    else:
        nband = max(np.ceil(nband*1.1), nband+4)

    # Increase number of bands based on the starting magnetization
    if nsppol==2 and spinat is not None:
        nband += np.ceil(max(np.sum(spinat, axis=0))/2.)

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

    if spin_mode=="polarized":
        multi[0].set_autospinat()

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, multi.pseudos, scf_electrons, multi[0].get('spinat', None))

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

    if spin_mode=="polarized":
        spinat_dict = multi[0].set_autospinat()
        multi[1].set_vars(spinat_dict)

    if electrons.nband is None:
        electrons.nband = _find_scf_nband(structure, multi.pseudos, electrons, multi[0].get('spinat', None))

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


def phonons_from_gsinput(gs_inp, ph_ngqpt=None, with_ddk=True, with_dde=True, with_bec=False, ph_tol=None, ddk_tol=None,
                         dde_tol=None):
    """
    Returns a :class:`AbinitInput` for performing phonon calculations.
    GS input + the input files for the phonon calculation.
    """

    gs_inp = gs_inp.deepcopy()

    gs_inp.pop_irdvars()

    if with_dde:
        with_ddk = True

    if with_bec:
        with_ddk = True
        with_dde = False

    multi = []

    if ph_ngqpt is None:
        ph_ngqpt = np.array(gs_inp["ngkpt"])

    qpoints = gs_inp.abiget_ibz(ngkpt=ph_ngqpt, shiftk=(0,0,0), kptopt=1).points


    # Build the input files for the q-points in the IBZ.
    # Response-function calculation for phonons.
    for qpt in qpoints:
        if np.allclose(qpt, 0):
            if with_ddk:
                multi_ddk = gs_inp.make_ddk_inputs(ddk_tol)
                multi_ddk.add_tags(DDK)
                multi.extend(multi_ddk)
            if with_dde:
                multi_dde = gs_inp.make_dde_inputs(dde_tol)
                multi_dde.add_tags(DDE)
                multi.extend(multi_dde)
            elif with_bec:
                multi_bec = gs_inp.make_bec_inputs(ph_tol)
                multi_bec.add_tags(BEC)
                multi.extend(multi_bec)
                continue

        multi_ph_q = gs_inp.make_ph_inputs_qpoint(qpt, ph_tol)
        multi_ph_q.add_tags(PH_Q_PERT)
        multi.extend(multi_ph_q)

    multi = MultiDataset.from_inputs(multi)
    multi.add_tags(PHONON)

    #FIXME for the time being there could be problems in mergeddb if the kpoints grid is gamma centered or if
    # if the grid is odd. Remove when mergeddb is fixed
    multi.set_vars(kptopt=3)

    return multi


def scf_piezo_elastic_inputs(structure, pseudos, kppa, ecut=None, pawecutdg=None, scf_nband=None,
                             accuracy="normal", spin_mode="polarized",
                             smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, ddk_tol=None, rf_tol=None):

    """
    Returns a :class:`AbinitInput` for performing elastic and piezoelectric constants calculations.
    GS input + the input files for the elastic and piezoelectric constants calculation.

    Args:
        structure: :class:`Structure` object.
        pseudos: List of filenames or list of :class:`Pseudo` objects or :class:`PseudoTable` object.
        kppa: Defines the sampling used for the SCF run.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the
            pseudos according to accuracy)
        scf_nband: Number of bands for SCF run. If scf_nband is None, nband is automatically initialized
            from the list of pseudos, the structure and the smearing option.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
        ddk_tol
    """
    # Build the input file for the GS run.
    gs_inp = AbinitInput(structure=structure, pseudos=pseudos)

    # Set the cutoff energies.
    gs_inp.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, gs_inp.pseudos))

    ksampling = aobj.KSampling.automatic_density(gs_inp.structure, kppa, chksymbreak=0, shifts=(0.0, 0.0, 0.0))
    gs_inp.set_vars(ksampling.to_abivars())
    gs_inp.set_vars(tolvrs=1.0e-18)

    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm,
                                   charge=charge, nband=None, fband=None)

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, gs_inp.pseudos, scf_electrons)
    gs_inp.set_vars(scf_electrons.to_abivars())

    all_inps = [gs_inp]

    # Add the ddk input
    ddk_inp = gs_inp.deepcopy()

    ddk_inp.set_vars(
                rfelfd=2,             # Activate the calculation of the d/dk perturbation
                rfdir=(1,1,1),        # All directions
                nqpt=1,               # One wavevector is to be considered
                qpt=(0, 0, 0),        # q-wavevector.
                kptopt=2,             # Take into account time-reversal symmetry.
                iscf=-3,              # The d/dk perturbation must be treated in a non-self-consistent way
            )
    if ddk_tol is None:
        ddk_tol = {"tolwfr": 1.0e-20}

    if len(ddk_tol) != 1 or any(k not in _tolerances for k in ddk_tol):
        raise ValueError("Invalid tolerance: {}".format(ddk_tol))
    ddk_inp.pop_tolerances()
    ddk_inp.set_vars(ddk_tol)

    ddk_inp.add_tags(DDK)
    all_inps.append(ddk_inp)

    # Add the Response Function calculation
    rf_inp = gs_inp.deepcopy()

    rf_inp.set_vars(rfphon=1,                          # Atomic displacement perturbation
                    rfatpol=(1,len(gs_inp.structure)), # Perturbation of all atoms
                    rfstrs=3,                          # Do the strain perturbations
                    rfdir=(1,1,1),                     # All directions
                    nqpt=1,                            # One wavevector is to be considered
                    qpt=(0, 0, 0),                     # q-wavevector.
                    kptopt=2,                          # Take into account time-reversal symmetry.
                    iscf=7,                            # The d/dk perturbation must be treated in a non-self-consistent way
                    )

    if rf_tol is None:
        rf_tol = {"tolvrs": 1.0e-12}

    if len(rf_tol) != 1 or any(k not in _tolerances for k in rf_tol):
        raise ValueError("Invalid tolerance: {}".format(rf_tol))
    rf_inp.pop_tolerances()
    rf_inp.set_vars(rf_tol)

    rf_inp.add_tags([DFPT, STRAIN])
    all_inps.append(rf_inp)

    return MultiDataset.from_inputs(all_inps)


def scf_input(structure, pseudos, kppa=None, ecut=None, pawecutdg=None, nband=None, accuracy="normal",
              spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
              shift_mode="Monkhorst-Pack"):
    structure = Structure.as_structure(structure)

    abinit_input = AbinitInput(structure, pseudos)

    # Set the cutoff energies.
    abinit_input.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, abinit_input.pseudos))

    # SCF calculation.
    kppa = _DEFAULTS.get("kppa") if kppa is None else kppa
    shifts = (0.5, 0.5, 0.5) if shift_mode[0].lower() == "m" else (0.0, 0.0, 0.0)
    scf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=shifts)
    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm,
                                   charge=charge, nband=nband, fband=None)

    if spin_mode=="polarized":
        abinit_input.set_autospinat()

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, abinit_input.pseudos, scf_electrons,abinit_input.get('spinat', None))

    abinit_input.set_vars(scf_ksampling.to_abivars())
    abinit_input.set_vars(scf_electrons.to_abivars())
    abinit_input.set_vars(_stopping_criterion("scf", accuracy))

    return abinit_input


def ebands_from_gsinput(gsinput, nband=None, ndivsm=15, accuracy="normal"):
    """

    :param gsinput:
    :param nband:
    :param ndivsm:
    :param accuracy:
    :return: AbinitInput
    """
    # create a copy to avoid messing with the previous input
    bands_input = gsinput.deepcopy()

    bands_input.pop_irdvars()

    nscf_ksampling = aobj.KSampling.path_from_structure(ndivsm, gsinput.structure)
    if nband is None:
        nband = gsinput.get("nband", gsinput.structure.num_valence_electrons(gsinput.pseudos)) + 10

    bands_input.set_vars(nscf_ksampling.to_abivars())
    bands_input.set_vars(nband=nband, iscf=-2)
    bands_input.set_vars(_stopping_criterion("nscf", accuracy))

    return bands_input


def ioncell_relax_from_gsinput(gsinput, accuracy="normal"):

    ioncell_input = gsinput.deepcopy()

    ioncell_input.pop_irdvars()

    ioncell_relax = aobj.RelaxationMethod.atoms_and_cell(atoms_constraints=None)
    ioncell_input.set_vars(ioncell_relax.to_abivars())
    ioncell_input.set_vars(_stopping_criterion("relax", accuracy))

    return ioncell_input


def hybrid_oneshot_input(gsinput, functional="hse06", ecutsigx=None, gw_qprange=1):

    hybrid_input = gsinput.deepcopy()

    hybrid_input.pop_irdvars()

    functional = functional.lower()
    if functional == 'hse06':
        gwcalctyp = 115
        icutcoul = 5
        rcut = 9.090909
    elif functional == 'pbe0':
        gwcalctyp = 215
        icutcoul = 6
        rcut = 0.
    elif functional == 'b3lyp':
        gwcalctyp = 315
        icutcoul = 6
        rcut = 0.
    else:
        raise ValueError("Unknow functional {}.".format(functional))

    ecut = hybrid_input['ecut']
    ecutsigx = ecutsigx or 2*ecut

    hybrid_input.set_vars(optdriver=4, gwcalctyp=gwcalctyp, gw_nstep=1, gwpara=2, icutcoul=icutcoul, rcut=rcut,
                          gw_qprange=gw_qprange, ecutwfn=ecut*0.995, ecutsigx=ecutsigx)

    return hybrid_input


def scf_for_phonons(structure, pseudos, kppa=None, ecut=None, pawecutdg=None, nband=None, accuracy="normal",
                    spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                    shift_mode="Symmetric"):
    abiinput = scf_input(structure=structure, pseudos=pseudos, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg, nband=nband,
                         accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
                         scf_algorithm=scf_algorithm, shift_mode=shift_mode)
    # set symmetrized k-point
    if shift_mode[0].lower() == 's':
        # need to convert to abipy structure to get the calc_shiftk method
        structure = Structure.from_sites(structure)
        shiftk = structure.calc_shiftk()
        abiinput.set_vars(shiftk=shiftk, nshiftk=len(shiftk))

    # enforce symmetries and add a buffer of bands to ease convergence with tolwfr
    abiinput.set_vars(chksymbreak=1, nbdbuf=4, tolwfr=1.e-22)

    return abiinput

#FIXME if the pseudos are passed as a PseudoTable the whole table will be serialized,
# it would be better to filter on the structure elements
class InputFactory(object):
    factory_function = None
    input_required = True

    def __init__(self, *args, **kwargs):
        if self.factory_function is None:
            raise NotImplementedError('The factory function should be specified')

        self.args = args
        self.kwargs = kwargs

    def build_input(self, previous_input=None):
        # make a copy to pop additional parameteres
        kwargs = dict(self.kwargs)
        decorators = kwargs.pop('decorators', [])
        if not isinstance(decorators, (list, tuple)):
            decorators = [decorators]
        extra_abivars = kwargs.pop('extra_abivars', {})
        if self.input_required:
            if not previous_input:
                raise ValueError('An input is required for factory function {}.'.format(self.factory_function.__name__))
            abiinput = self.factory_function(previous_input, *self.args, **kwargs)
        else:
            abiinput = self.factory_function(*self.args, **kwargs)

        for d in decorators:
            abiinput = d(abiinput)
        abiinput.set_vars(extra_abivars)

        return abiinput

    @pmg_serialize
    def as_dict(self):
        # sanitize to avoid numpy arrays and serialize PMGSonable objects
        return jsanitize(dict(args=self.args, kwargs=self.kwargs), strict=True)

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        return cls(*dec.process_decoded(d['args']), **dec.process_decoded(d['kwargs']))


class BandsFromGsFactory(InputFactory):
    factory_function = staticmethod(ebands_from_gsinput)


class IoncellRelaxFromGsFactory(InputFactory):
    factory_function = staticmethod(ioncell_relax_from_gsinput)


class HybridOneShotFromGsFactory(InputFactory):
    factory_function = staticmethod(hybrid_oneshot_input)


class ScfFactory(InputFactory):
    factory_function = staticmethod(scf_input)
    input_required = False

class ScfForPhononsFactory(InputFactory):
    factory_function = staticmethod(scf_for_phonons)
    input_required = False


class PhononsFromGsFactory(InputFactory):
    factory_function = staticmethod(phonons_from_gsinput)
