# coding: utf-8
"""Factory functions for Abinit input files """
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import pymatgen.io.abinit.abiobjects as aobj

from enum import Enum
from collections import namedtuple
from monty.collections import AttrDict
from monty.string import is_string
from monty.json import jsanitize, MontyDecoder
try:
    from pymatgen.util.serialization import pmg_serialize
except ImportError:
    from pymatgen.serializers.json_coders import pmg_serialize
from abipy.flowtk import PseudoTable
from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput, MultiDataset
from abipy.abio.input_tags import *
from monty.json import MSONable

import logging
logger = logging.getLogger(__file__)


__all__ = [
    "gs_input",
    "ebands_input",
    "phonons_from_gsinput",
    "g0w0_with_ppmodel_inputs",
    "g0w0_convergence_inputs",
    "bse_with_mdf_inputs",
    "ion_ioncell_relax_input",
    "ion_ioncell_relax_and_ebands_input",
    "scf_phonons_inputs",
    "piezo_elastic_inputs_from_gsinput",
    "scf_piezo_elastic_inputs",
    "scf_for_phonons",
    "dte_from_gsinput",
    "dfpt_from_gsinput",
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

class ShiftMode(Enum):
    """
    Class defining the mode to be used for the shifts.
    G: Gamma centered
    M: Monkhorst-Pack ((0.5, 0.5, 0.5))
    S: Symmetric. Respects the chksymbreak with multiple shifts
    O: OneSymmetric. Respects the chksymbreak with a single shift (as in 'S' if a single shift is given, gamma
        centered otherwise.
    """
    GammaCentered = 'G'
    MonkhorstPack = 'M'
    Symmetric = 'S'
    OneSymmetric = 'O'

    @classmethod
    def from_object(cls, obj):
        """
        Returns an instance of ShiftMode based on the type of object passed. Converts strings to ShiftMode depending
        on the iniital letter of the string. G for GammaCenterd, M for MonkhorstPack, S for Symmetric, O for OneSymmetric.
        Case insensitive.
        """
        if isinstance(obj, cls):
            return obj
        elif is_string(obj):
            return cls(obj[0].upper())
        else:
            raise TypeError('The object provided is not handled: type' % type(obj))


def _stopping_criterion(runlevel, accuracy):
    """Return the stopping criterion for this runlevel with the given accuracy."""
    tolname = _runl2tolname[runlevel]
    return {tolname: getattr(_tolerances[tolname], accuracy)}


def _find_ecut_pawecutdg(ecut, pawecutdg, pseudos, accuracy):
    """Return a |AttrDict| with the value of ``ecut`` and ``pawecutdg``."""
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
    """Find the value of ``nband``."""
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
    if nsppol == 2 and spinat is not None:
        nband += np.ceil(max(np.sum(spinat, axis=0))/2.)

    # Force even nband (easier to divide among procs, mandatory if nspinor == 2)
    nband += nband % 2
    return int(nband)


def _get_shifts(shift_mode, structure):
    """
    Gives the shifts based on the selected shift mode and on the symmetry of the structure.
    G: Gamma centered
    M: Monkhorst-Pack ((0.5, 0.5, 0.5))
    S: Symmetric. Respects the chksymbreak with multiple shifts
    O: OneSymmetric. Respects the chksymbreak with a single shift (as in 'S' if a single shift is given, gamma
        centered otherwise.

    Note: for some cases (e.g. body centered tetragonal), both the Symmetric and OneSymmetric may fail to satisfy the
        ``chksymbreak`` condition (Abinit input variable).
    """
    if shift_mode == ShiftMode.GammaCentered:
        return ((0, 0, 0))
    elif shift_mode == ShiftMode.MonkhorstPack:
        return ((0.5, 0.5, 0.5))
    elif shift_mode == ShiftMode.Symmetric:
        structure = Structure.from_sites(structure)
        return structure.calc_shiftk()
    elif shift_mode == ShiftMode.OneSymmetric:
        structure = Structure.from_sites(structure)
        shifts = structure.calc_shiftk()
        if len(shifts) == 1:
            return shifts
        else:
            return ((0, 0, 0))
    else:
        raise ValueError("shift_mode `%s `not valid." % str(shift_mode))


def gs_input(structure, pseudos,
             kppa=None, ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
             smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
    """
    Returns a |AbinitInput| for ground-state calculation.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
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
                 kppa=kppa, ndivsm=0,
                 ecut=ecut, pawecutdg=pawecutdg, scf_nband=scf_nband, accuracy=accuracy, spin_mode=spin_mode,
                 smearing=smearing, charge=charge, scf_algorithm=scf_algorithm)

    return multi[0]


def ebands_input(structure, pseudos,
                 kppa=None, nscf_nband=None, ndivsm=15,
                 ecut=None, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode="polarized",
                 smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, dos_kppa=None):
    """
    Returns a |MultiDataset| object for band structure calculations.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
        kppa: Defines the sampling used for the SCF run. Defaults to 1000 if not given.
        nscf_nband: Number of bands included in the NSCF run. Set to scf_nband + 10 if None.
        ndivsm: Number of divisions used to sample the smallest segment of the k-path.
            if 0, only the GS input is returned in multi[0].
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
    multi.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos, accuracy))

    # SCF calculation.
    kppa = _DEFAULTS.get("kppa") if kppa is None else kppa
    scf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm,
                                   charge=charge, nband=scf_nband, fband=None)

    if spin_mode == "polarized":
        multi[0].set_autospinat()

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, multi.pseudos, scf_electrons, multi[0].get('spinat', None))

    multi[0].set_vars(scf_ksampling.to_abivars())
    multi[0].set_vars(scf_electrons.to_abivars())
    multi[0].set_vars(_stopping_criterion("scf", accuracy))
    if ndivsm == 0: return multi

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
                            smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, shift_mode='Monkhorst-pack'):
    """
    Returns a |MultiDataset| for a structural relaxation. The first dataset optmizes the
    atomic positions at fixed unit cell. The second datasets optimizes both ions and unit cell parameters.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
        kppa: Defines the sampling used for the Brillouin zone.
        nband: Number of bands included in the SCF run.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for the solution of the SCF cycle.
    """
    structure = Structure.as_structure(structure)
    multi = MultiDataset(structure, pseudos, ndtset=2)

    # Set the cutoff energies.
    multi.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos, accuracy))

    kppa = _DEFAULTS.get("kppa") if kppa is None else kppa

    shift_mode = ShiftMode.from_object(shift_mode)
    shifts = _get_shifts(shift_mode, structure)
    ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=shifts)
    electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm,
                               charge=charge, nband=nband, fband=None)

    if spin_mode == "polarized":
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
    Returns a |MultiDataset| for a structural relaxation followed by a band structure run.
    The first dataset optimizes the atomic positions at fixed unit cell.
    The second datasets optimizes both ions and unit cell parameters.
    The other datasets perform a band structure calculation.

    .. warning::

        Client code is responsible for propagating the relaxed structure obtained with the
        second dataset to the inputs used for the band structure calculation.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
        kppa: Defines the sampling used for the Brillouin zone.
        nband: Number of bands included in the SCF run.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.

    Returns: |MultiDataset| object
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
                             ecut=None, pawecutdg=None, shifts=(0.0, 0.0, 0.0),
                             accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
                             ppmodel="godby", charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None,
                             sigma_nband=None, gw_qprange=1):
    """
    Returns a |MultiDataset| object that performs G0W0 calculations with the plasmon pole approximation.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
        kppa: Defines the sampling used for the SCF run.
        nscf_nband: Number of bands included in the NSCF run.
        ecuteps: Cutoff energy [Ha] for the screening matrix.
        ecutsigx: Cutoff energy [Ha] for the exchange part of the self-energy.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized
            from the pseudos according to accuracy)
        shifts: Shifts for k-mesh.
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

    .. versionchanged: 0.3

        The default value of ``shifts`` changed in v0.3 from (0.5, 0.5, 0.5) to (0.0, 0.0, 0.0).
    """

    structure = Structure.as_structure(structure)
    multi = MultiDataset(structure, pseudos, ndtset=4)

    # Set the cutoff energies.
    multi.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos, accuracy))

    scf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=shifts)
    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm,
                                   charge=charge, nband=None, fband=None)

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, multi.pseudos, scf_electrons)

    multi[0].set_vars(scf_ksampling.to_abivars())
    multi[0].set_vars(scf_electrons.to_abivars())
    multi[0].set_vars(_stopping_criterion("scf", accuracy))

    nscf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=shifts)
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

    # Sigma.
    if sigma_nband is None: sigma_nband = nscf_nband
    self_energy = aobj.SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening,
                                  gw_qprange=gw_qprange, ppmodel=ppmodel)

    multi[3].set_vars(nscf_ksampling.to_abivars())
    multi[3].set_vars(nscf_electrons.to_abivars())
    multi[3].set_vars(self_energy.to_abivars())
    multi[3].set_vars(_stopping_criterion("sigma", accuracy)) # Dummy

    # TODO: Cannot use istwfk != 1.
    multi.set_vars(istwfk="*1")

    return multi


def g0w0_convergence_inputs(structure, pseudos, kppa, nscf_nband, ecuteps, ecutsigx, scf_nband, ecut,
                            accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
                            response_models=None, charge=0.0, scf_algorithm=None, inclvkb=2,
                            gw_qprange=1, gamma=True, nksmall=None, extra_abivars=None):
    """
    Returns a |MultiDataset| object to generate a G0W0 work for the given the material.
    See also :cite:`Setten2017`.

    Args:
        structure: |Structure| object
        pseudos: List of |Pseudo| objects.
        kppa: k points per reciprocal atom.
        scf_nband: number of scf bands
        ecut: ecut for all calcs that that are not ecut convergence  cals at scf level
        scf_ Defines the sampling used for the SCF run.
        nscf_nband: a list of number of bands included in the screening and sigmaruns.
            The NSCF run will be done with the maximum.
        ecuteps: list of Cutoff energy [Ha] for the screening matrix.
        ecutsigx: Cutoff energy [Ha] for the exchange part of the self-energy.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
        inclvkb: Treatment of the dipole matrix elements (see abinit variable).
        response_models: List of response models
        gw_qprange: selectpr for the qpoint mesh
        gamma: is true a gamma centered mesh is enforced
        nksmall: Kpoint division for additional band and dos calculations
        extra_abivars: Dictionary with extra variables passed to ABINIT for all tasks.

    extra abivars that are provided with _s appended will be take as a list of values to be tested a scf level
    """
    if extra_abivars is None:
        extra_abivars = {}

    if response_models is None:
        response_models = ["godby"]

    scf_diffs = []

    keys = list(extra_abivars.keys())
    #for k in extra_abivars.keys():
    for k in keys:
        if k[-2:] == '_s':
            var = k[:len(k)-2]
            values = extra_abivars.pop(k)
            # to_add.update({k: values[-1]})
            for value in values:
                diff_abivars = dict()
                diff_abivars[var] = value
                if pseudos.allpaw and var == 'ecut':
                    diff_abivars['pawecutdg'] = diff_abivars['ecut'] * 2
                scf_diffs.append(diff_abivars)

    extra_abivars_all = dict(
        ecut=ecut,
        paral_kgb=1,
        istwfk="*1",
        timopt=-1,
        nbdbuf=8,
    )

    extra_abivars_all.update(extra_abivars)

    if pseudos.allpaw:
        extra_abivars_all['pawecutdg'] = extra_abivars_all['ecut'] * 2

    extra_abivars_gw = dict(
        inclvkb=2,
        symsigma=1,
        gwpara=2,
        gwmem='10',
        prtsuscep=0
    )

    # all these too many options are for development only the current idea for the final version is
    #if gamma:
    #    scf_ksampling = aobj.KSampling.automatic_density(structure=structure, kppa=10000, chksymbreak=0, shifts=(0, 0, 0))
    #    nscf_ksampling = aobj.KSampling.gamma_centered(kpts=(2, 2, 2))
    #    if kppa <= 13:
    #        nscf_ksampling = aobj.KSampling.gamma_centered(kpts=(scf_kppa, scf_kppa, scf_kppa))
    #    else:
    #        nscf_ksampling = aobj.KSampling.automatic_density(structure, scf_kppa, chksymbreak=0, shifts=(0, 0, 0))
    #else:
    #    scf_ksampling = aobj.KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
    #    nscf_ksampling = aobj.KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)

    if gamma:
        if kppa == 1:
            scf_ksampling = aobj.KSampling.gamma_centered(kpts=(1, 1, 1))
            nscf_ksampling = aobj.KSampling.gamma_centered(kpts=(1, 1, 1))
        elif kppa == 2:
            scf_ksampling = aobj.KSampling.gamma_centered(kpts=(2, 2, 2))
            nscf_ksampling = aobj.KSampling.gamma_centered(kpts=(2, 2, 2))
        elif kppa < 0:
            scf_ksampling = aobj.KSampling.gamma_centered(kpts=(-kppa, -kppa, -kppa))
            nscf_ksampling = aobj.KSampling.gamma_centered(kpts=(2, 2, 2))
        elif kppa <= 13:
            scf_ksampling = aobj.KSampling.gamma_centered(kpts=(kppa, kppa, kppa))
            nscf_ksampling = aobj.KSampling.gamma_centered(kpts=(kppa, kppa, kppa))
        else:
            scf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=(0, 0, 0))
            nscf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=(0, 0, 0))
    else:
        # this is the original behaviour before the development of the gwwrapper
        scf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
        nscf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)

    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm,
                                   charge=charge, nband=scf_nband, fband=None)
    nscf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
                                    charge=charge, nband=max(nscf_nband), fband=None)

    multi_scf = MultiDataset(structure, pseudos, ndtset=max(1, len(scf_diffs)))

    multi_scf.set_vars(scf_ksampling.to_abivars())
    multi_scf.set_vars(scf_electrons.to_abivars())
    multi_scf.set_vars(extra_abivars_all)
    multi_scf.set_vars(_stopping_criterion(runlevel="scf", accuracy=accuracy))
    multi_scf.set_vars(extra_abivars)

    for variables, abinput in zip(scf_diffs, multi_scf):
        abinput.set_vars(variables)

    scf_inputs = multi_scf.split_datasets()

    # create nscf inputs
    ndtset = 3 if nksmall is not None else 1
    nscf_multi = MultiDataset(structure=structure, pseudos=pseudos, ndtset=ndtset)

    nscf_multi.set_vars(nscf_electrons.to_abivars())
    nscf_multi.set_vars(extra_abivars_all)
    nscf_multi.set_vars(_stopping_criterion(runlevel="nscf", accuracy=accuracy))

    nscf_multi[-1].set_vars(nscf_ksampling.to_abivars())

    if nksmall is not None:
        # if nksmall add bandstructure and dos calculations as well
        logger.info('added band structure calculation')
        bands_ksampling = aobj.KSampling.path_from_structure(ndivsm=nksmall, structure=structure)
        dos_ksampling = aobj.KSampling.automatic_density(structure=structure, kppa=2000)
        nscf_multi[0].set_vars(bands_ksampling.to_abivars())
        nscf_multi[0].set_vars({'chksymbreak': 0})
        nscf_multi[1].set_vars(dos_ksampling.to_abivars())
        nscf_multi[1].set_vars({'chksymbreak': 0})

    nscf_inputs = nscf_multi.split_datasets()

    # create screening and sigma inputs

    #if scr_nband is None:
    #   scr_nband = nscf_nband_nscf
    #if sigma_nband is None:
    #     sigma_nband = nscf_nband_nscf

    if 'cd' in response_models:
        hilbert = aobj.HilbertTransform(nomegasf=100, domegasf=None, spmeth=1, nfreqre=None, freqremax=None, nfreqim=None,
                                        freqremin=None)
    scr_inputs = []
    sigma_inputs = []
    #print("ecuteps", ecuteps, "nscf_nband", nscf_nband)

    for response_model in response_models:
        for ecuteps_v in ecuteps:
            for nscf_nband_v in nscf_nband:
                scr_nband = nscf_nband_v
                sigma_nband = nscf_nband_v
                multi = MultiDataset(structure, pseudos, ndtset=2)
                multi.set_vars(nscf_ksampling.to_abivars())
                multi.set_vars(nscf_electrons.to_abivars())
                multi.set_vars(extra_abivars_all)
                multi.set_vars(extra_abivars_gw)
                if response_model == 'cd':
                    screening = aobj.Screening(ecuteps_v, scr_nband, w_type="RPA", sc_mode="one_shot", hilbert=hilbert,
                                               ecutwfn=None, inclvkb=inclvkb)
                    self_energy = aobj.SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening)
                else:
                    ppmodel = response_model
                    screening = aobj.Screening(ecuteps_v, scr_nband, w_type="RPA", sc_mode="one_shot",
                                               hilbert=None, ecutwfn=None, inclvkb=inclvkb)
                    self_energy = aobj.SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening,
                                                  gw_qprange=gw_qprange, ppmodel=ppmodel)
                multi[0].set_vars(screening.to_abivars())
                multi[0].set_vars(_stopping_criterion("screening", accuracy))  # Dummy
                multi[1].set_vars(self_energy.to_abivars())
                multi[1].set_vars(_stopping_criterion("sigma", accuracy))  # Dummy

                scr_input, sigma_input = multi.split_datasets()
                scr_inputs.append(scr_input)
                sigma_inputs.append(sigma_input)

    return scf_inputs, nscf_inputs, scr_inputs, sigma_inputs


def bse_with_mdf_inputs(structure, pseudos,
                        scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
                        ecuteps, bs_loband, bs_nband, mbpt_sciss, mdf_epsinf,
                        ecut=None, pawecutdg=None,
                        exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="polarized",
                        smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None):
    """
    Returns a |MultiDataset| object that performs a GS + NSCF + Bethe-Salpeter calculation.
    The self-energy corrections are approximated with the scissors operator.
    The screening is modeled with the model dielectric function.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
        scf_kppa: Defines the sampling used for the SCF run.
        nscf_nband: Number of bands included in the NSCF run.
        nscf_ngkpt: Divisions of the k-mesh used for the NSCF and the BSE run.
        nscf_shiftk: Shifts used for the NSCF and the BSE run.
        ecuteps: Cutoff energy [Ha] for the screening matrix.
        bs_loband: Index of the first occupied band included the e-h basis set
            (ABINIT convention i.e. first band starts at 1).
            Can be scalar or array of shape (nsppol,)
        bs_nband: Highest band idex used for the construction of the e-h basis set.
        mbpt_sciss: Scissor energy in Hartree.
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
    d = _find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos, accuracy)
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
    exc_ham = aobj.ExcHamiltonian(bs_loband, bs_nband, mbpt_sciss, coulomb_mode="model_df", ecuteps=ecuteps,
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
    # TODO: Please check the unused variables in the function
    """
    Returns a list of input files for performing phonon calculations.
    GS input + the input files for the phonon calculation.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
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
    gs_inp.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, gs_inp.pseudos, accuracy))

    ksampling = aobj.KSampling.automatic_density(gs_inp.structure, kppa, chksymbreak=0)
    gs_inp.set_vars(ksampling.to_abivars())
    gs_inp.set_vars(tolvrs=1.0e-18)

    # Get the qpoints in the IBZ. Note that here we use a q-mesh with ngkpt=(4,4,4) and shiftk=(0,0,0)
    # i.e. the same parameters used for the k-mesh in gs_inp.
    qpoints = gs_inp.abiget_ibz(ngkpt=(4, 4, 4), shiftk=(0, 0, 0), kptopt=1).points
    #print("get_ibz qpoints:", qpoints)

    # Build the input files for the q-points in the IBZ.
    #ph_inputs = MultiDataset(gs_inp.structure, pseudos=gs_inp.pseudos, ndtset=len(qpoints))

    ph_inputs = MultiDataset.replicate_input(gs_inp, ndtset=len(qpoints))

    for ph_inp, qpt in zip(ph_inputs, qpoints):
        # Response-function calculation for phonons.
        ph_inp.set_vars(
            rfphon=1,        # Will consider phonon-type perturbation
            nqpt=1,          # One wavevector is to be considered
            qpt=qpt,         # This wavevector is q=0 (Gamma)
            tolwfr=1.0e-20,
            kptopt=3,        # TODO: One could use symmetries for Gamma.
        )
            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            #kptopt   2      # Automatic generation of k points, taking

        irred_perts = ph_inp.abiget_irred_phperts()
        # TODO irred_perts is not used ??

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


def phonons_from_gsinput(gs_inp, ph_ngqpt=None, qpoints=None, with_ddk=True, with_dde=True, with_bec=False,
                         ph_tol=None, ddk_tol=None, dde_tol=None, wfq_tol=None, qpoints_to_skip=None, manager=None):
    """
    Returns a list of inputs in the form of a MultiDataset to perform phonon calculations, based on
    a ground state |AbinitInput|.
    It will determine if WFQ files should be calculated for some q points and add the NSCF AbinitInputs to the set.
    The inputs have the following tags, according to their function: "ddk", "dde", "nscf", "ph_q_pert".
    All of them have the tag "phonon".

    Args:
        gs_inp: an |AbinitInput| representing a ground state calculation, likely the SCF performed to get the WFK.
        ph_ngqpt: a list of three integers representing the gamma centered q-point grid used for the calculation.
            If None and qpoint==None the ngkpt value present in the gs_input will be used.
            Incompatible with qpoints.
        qpoints: a list of coordinates of q points in reduced coordinates for which the phonon perturbations will
            be calculated. Incompatible with ph_ngqpt.
        with_ddk: If True, if Gamma is included in the list of qpoints it will add inputs for the calculations of
            the DDK.
        with_dde: If True, if Gamma is included in the list of qpoints it will add inputs for the calculations of
            the DDE. Automatically sets with_ddk=True.
        with_bec: If Truem if Gamma is included in the list of qpoints the DDE will be calculated in the same
            input as the phonons. This will allow to determine the BECs.
            Automatically sets with_ddk=True and with_dde=False.
        ph_tol: a dictionary with a single key defining the type of tolerance used for the phonon calculations and
            its value. Default: {"tolvrs": 1.0e-10}.
        ddk_tol: a dictionary with a single key defining the type of tolerance used for the DDK calculations and
            its value. Default: {"tolwfr": 1.0e-22}.
        dde_tol: a dictionary with a single key defining the type of tolerance used for the DDE calculations and
            its value. Default: {"tolvrs": 1.0e-10}.
        wfq_tol: a dictionary with a single key defining the type of tolerance used for the NSCF calculations of
            the WFQ and its value. Default {"tolwfr": 1.0e-22}.
        qpoints_to_skip: a list of coordinates of q points in reduced coordinates that will be skipped.
            Useful when calculating multiple grids for the same system to avoid duplicate calculations.
            If a DDB needs to be extended with more q points use e.g. ddb.qpoints.to_array().
        manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
    """
    gs_inp = gs_inp.deepcopy()
    gs_inp.pop_irdvars()

    if with_dde:
        with_ddk = True

    if with_bec:
        with_ddk = True
        with_dde = False

    if ph_tol is None:
        ph_tol = {"tolvrs": 1.0e-10}

    if ddk_tol is None:
        ddk_tol = {"tolwfr": 1.0e-22}

    if dde_tol is None:
        dde_tol = {"tolvrs": 1.0e-10}

    if wfq_tol is None:
        wfq_tol = {"tolwfr": 1.0e-22}

    multi = []

    if qpoints is not None and ph_ngqpt is not None:
        raise ValueError("ph_ngqpt and qpoints can't be used together")

    if qpoints is None:
        if ph_ngqpt is None:
            ph_ngqpt = np.array(gs_inp["ngkpt"])
        else:
            ph_ngqpt = np.array(ph_ngqpt)

        qpoints = gs_inp.abiget_ibz(ngkpt=ph_ngqpt, shiftk=(0, 0, 0), kptopt=1, manager=manager).points

    if qpoints_to_skip:
        preserved_qpoints = []
        for q in qpoints:
            if not any(np.allclose(q, ddb_q) for ddb_q in qpoints_to_skip):
                preserved_qpoints.append(q)
        qpoints = np.array(preserved_qpoints)

    if ph_ngqpt is None or any(gs_inp["ngkpt"] % ph_ngqpt != 0):
        # find which q points are needed and build nscf inputs to calculate the WFQ
        kpts = gs_inp.abiget_ibz(shiftk=(0, 0, 0), kptopt=3, manager=manager).points.tolist()
        nscf_qpt = []
        for q in qpoints:
            if list(q) not in kpts:
                nscf_qpt.append(q)
        if nscf_qpt:
            multi_nscf = MultiDataset.replicate_input(gs_inp, len(nscf_qpt))
            multi_nscf.set_vars(kptopt=3, nqpt=1, iscf=-2)
            if wfq_tol:
                multi_nscf.set_vars(**wfq_tol)
            else:
                multi_nscf.set_vars(tolwfr=1e-22)
            for q, nscf_inp in zip(nscf_qpt, multi_nscf):
                nscf_inp.set_vars(qpt=q)

            multi_nscf.add_tags(NSCF)

            multi.extend(multi_nscf)

    # Build the input files for the q-points in the IBZ.
    # Response-function calculation for phonons.
    for qpt in qpoints:
        if np.allclose(qpt, 0):
            if with_ddk:
                multi_ddk = gs_inp.make_ddk_inputs(tolerance=ddk_tol)
                multi_ddk.add_tags(DDK)
                multi.extend(multi_ddk)
            if with_dde:
                multi_dde = gs_inp.make_dde_inputs(dde_tol, manager=manager)
                multi_dde.add_tags(DDE)
                multi.extend(multi_dde)
            elif with_bec:
                multi_bec = gs_inp.make_bec_inputs(ph_tol, manager=manager)
                multi_bec.add_tags(BEC)
                multi.extend(multi_bec)
                continue

        multi_ph_q = gs_inp.make_ph_inputs_qpoint(qpt, ph_tol)
        multi_ph_q.add_tags(PH_Q_PERT)
        multi.extend(multi_ph_q)

    multi = MultiDataset.from_inputs(multi)
    multi.add_tags(PHONON)

    return multi

def piezo_elastic_inputs_from_gsinput(gs_inp, ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False,
                                      manager=None):
    """
    Returns a |MultiDataset| for performing elastic and piezoelectric constants calculations.
    GS input + the input files for the elastic and piezoelectric constants calculation.

    Args:
        gs_inp: Ground State input to build piezo elastic inputs from.
        ddk_tol: Tolerance for the DDK calculation (i.e. {"tolwfr": 1.0e-20}).
        rf_tol: Tolerance for the Strain RF calculations (i.e. {"tolvrs": 1.0e-12}).
        ddk_split: Whether to split the DDK calculations.
        rf_split: whether to split the RF calculations.
        manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
    """
    # Ddk input(s)
    if ddk_split:
        multi = gs_inp.make_ddk_inputs(tolerance=ddk_tol)
    else:
        ddk_inp = gs_inp.deepcopy()

        ddk_inp.set_vars(
                    rfelfd=2,             # Activate the calculation of the d/dk perturbation
                    rfdir=(1,1,1),        # All directions
                    nqpt=1,               # One wavevector is to be considered
                    qpt=(0, 0, 0),        # q-wavevector.
                    kptopt=3,             # Take into account time-reversal symmetry.
                    iscf=-3,              # The d/dk perturbation must be treated in a non-self-consistent way
                    paral_kgb=0
                )
        if ddk_tol is None:
            ddk_tol = {"tolwfr": 1.0e-20}

        if len(ddk_tol) != 1 or any(k not in _tolerances for k in ddk_tol):
            raise ValueError("Invalid tolerance: {}".format(ddk_tol))
        ddk_inp.pop_tolerances()
        ddk_inp.set_vars(ddk_tol)
        # Adding buffer to help convergence ...
        if 'nbdbuf' not in ddk_inp:
            nbdbuf = max(int(0.1*ddk_inp['nband']), 4)
            ddk_inp.set_vars(nband=ddk_inp['nband']+nbdbuf, nbdbuf=nbdbuf)

        multi = MultiDataset.from_inputs([ddk_inp])
    multi.add_tags(DDK)

    # Response Function input(s)
    if rf_split:
        multi_rf = gs_inp.make_strain_perts_inputs(tolerance=rf_tol, manager=manager)
    else:
        rf_inp = gs_inp.deepcopy()

        rf_inp.set_vars(rfphon=1,                          # Atomic displacement perturbation
                        rfatpol=(1,len(gs_inp.structure)), # Perturbation of all atoms
                        rfstrs=3,                          # Do the strain perturbations
                        rfdir=(1,1,1),                     # All directions
                        nqpt=1,                            # One wavevector is to be considered
                        qpt=(0, 0, 0),                     # q-wavevector.
                        kptopt=3,                          # Take into account time-reversal symmetry.
                        iscf=7,                            # The rfstrs perturbation must be treated in a
                                                           #  self-consistent way
                        paral_kgb=0
                        )

        if rf_tol is None:
            rf_tol = {"tolvrs": 1.0e-12}

        if len(rf_tol) != 1 or any(k not in _tolerances for k in rf_tol):
            raise ValueError("Invalid tolerance: {}".format(rf_tol))
        rf_inp.pop_tolerances()
        rf_inp.set_vars(rf_tol)

        # Adding buffer to help convergence ...
        if 'nbdbuf' not in rf_inp:
            nbdbuf = max(int(0.1*rf_inp['nband']), 4)
            rf_inp.set_vars(nband=rf_inp['nband']+nbdbuf, nbdbuf=nbdbuf)

        multi_rf = MultiDataset.from_inputs([rf_inp])
    multi_rf.add_tags([DFPT, STRAIN])
    for inp in multi_rf:
        if inp.get('rfphon', 0) == 1:
            inp.add_tags(PHONON)

    multi.extend(multi_rf)

    return multi


def scf_piezo_elastic_inputs(structure, pseudos, kppa, ecut=None, pawecutdg=None, scf_nband=None,
                             accuracy="normal", spin_mode="polarized",
                             smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                             ddk_tol=None, rf_tol=None, ddk_split=False, rf_split=False):

    """
    Returns a |MultiDataset| for performing elastic and piezoelectric constants calculations.
    GS input + the input files for the elastic and piezoelectric constants calculation.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
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
        ddk_tol: Tolerance for the Ddk calculation (i.e. {"tolwfr": 1.0e-20}).
        rf_tol: Tolerance for the Strain RF calculations (i.e. {"tolvrs": 1.0e-12}).
        ddk_split: Whether to split the ddk calculations.
        rf_split: whether to split the RF calculations.
    """
    # Build the input file for the GS run.
    gs_inp = scf_input(structure=structure, pseudos=pseudos, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg,
                       nband=scf_nband, accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
                       scf_algorithm=scf_algorithm, shift_mode="Gamma-centered")

    # Adding buffer to help convergence ...
    nbdbuf = max(int(0.1*gs_inp['nband']), 4)
    gs_inp.set_vars(nband=gs_inp['nband']+nbdbuf, nbdbuf=nbdbuf)

    multi = MultiDataset.from_inputs([gs_inp])

    piezo_elastic_inputs = piezo_elastic_inputs_from_gsinput(gs_inp=gs_inp, ddk_tol=ddk_tol, rf_tol=rf_tol)

    multi.extend(piezo_elastic_inputs)

    return multi


def scf_input(structure, pseudos, kppa=None, ecut=None, pawecutdg=None, nband=None, accuracy="normal",
              spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
              shift_mode="Monkhorst-Pack"):
    """
    Returns an |AbinitInput| object for standard GS calculations.
    """
    structure = Structure.as_structure(structure)

    abinit_input = AbinitInput(structure, pseudos)

    # Set the cutoff energies.
    abinit_input.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, abinit_input.pseudos, accuracy))

    # SCF calculation.
    kppa = _DEFAULTS.get("kppa") if kppa is None else kppa
    shift_mode = ShiftMode.from_object(shift_mode)
    shifts = _get_shifts(shift_mode, structure)
    scf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=shifts)
    scf_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm,
                                   charge=charge, nband=nband, fband=None)

    if spin_mode == "polarized":
        abinit_input.set_autospinat()

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, abinit_input.pseudos, scf_electrons,
                                              abinit_input.get('spinat', None))

    abinit_input.set_vars(scf_ksampling.to_abivars())
    abinit_input.set_vars(scf_electrons.to_abivars())
    abinit_input.set_vars(_stopping_criterion("scf", accuracy))

    return abinit_input


def ebands_from_gsinput(gsinput, nband=None, ndivsm=15, accuracy="normal"):
    """
    Return an |AbinitInput| object to compute a band structure from a GS SCF input.

    Args:
        gsinput:
        nband:
        ndivsm:
        accuracy:

    Return: |AbinitInput|
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


def dos_from_gsinput(gsinput, dos_kppa, nband=None, accuracy="normal", pdos=False):

    # create a copy to avoid messing with the previous input
    dos_input = gsinput.deepcopy()
    dos_input.pop_irdvars()

    dos_ksampling = aobj.KSampling.automatic_density(dos_input.structure, dos_kppa, chksymbreak=0)
    dos_input.set_vars(dos_ksampling.to_abivars())
    dos_input.set_vars(iscf=-2, ionmov=0)
    dos_input.set_vars(_stopping_criterion("nscf", accuracy))

    if pdos:
        # FIXME
        raise NotImplementedError()

    return dos_input


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
        raise ValueError("Unknow functional {0}.".format(functional))

    ecut = hybrid_input['ecut']
    ecutsigx = ecutsigx or 2*ecut

    hybrid_input.set_vars(optdriver=4, gwcalctyp=gwcalctyp, gw_nstep=1, gwpara=2, icutcoul=icutcoul, rcut=rcut,
                          gw_qprange=gw_qprange, ecutwfn=ecut*0.995, ecutsigx=ecutsigx)

    return hybrid_input


def hybrid_scf_input(gsinput, functional="hse06", ecutsigx=None, gw_qprange=1):

    hybrid_input = hybrid_oneshot_input(gsinput=gsinput, functional=functional, ecutsigx=ecutsigx, gw_qprange=gw_qprange)

    hybrid_input['gwcalctyp'] += 10

    return hybrid_input


def scf_for_phonons(structure, pseudos, kppa=None, ecut=None, pawecutdg=None, nband=None, accuracy="normal",
                    spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                    shift_mode="Symmetric"):

    abiinput = scf_input(structure=structure, pseudos=pseudos, kppa=kppa, ecut=ecut, pawecutdg=pawecutdg, nband=nband,
                         accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, charge=charge,
                         scf_algorithm=scf_algorithm, shift_mode=shift_mode)

    nbdbuf = 4
    # with no smearing set the minimum number of bands plus some nbdbuf
    if smearing is None:
        nval = structure.num_valence_electrons(pseudos)
        nval -= abiinput['charge']
        nband = int(round(nval / 2) + nbdbuf)
        abiinput.set_vars(nband=nband)

    # enforce symmetries and add a buffer of bands to ease convergence with tolwfr
    abiinput.set_vars(chksymbreak=1, nbdbuf=nbdbuf, tolwfr=1.e-22)

    return abiinput


def dte_from_gsinput(gs_inp, use_phonons=True, ph_tol=None, ddk_tol=None, dde_tol=None,
                     skip_dte_permutations=False, manager=None):
    """
    Returns a list of inputs in the form of a |MultiDataset| to perform calculations of non-linear properties, based on
    a ground state AbinitInput.

    The inputs have the following tags, according to their function: "ddk", "dde", "ph_q_pert" and "dte".
    All of them have the tag "dfpt".

    Args:
        gs_inp: an |AbinitInput| representing a ground state calculation, likely the SCF performed to get the WFK.
        use_phonons: determine wether the phonon perturbations at gamma should be included or not
        ph_tol: a dictionary with a single key defining the type of tolerance used for the phonon calculations and
            its value. Default: {"tolvrs": 1.0e-22}.
        ddk_tol: a dictionary with a single key defining the type of tolerance used for the DDK calculations and
            its value. Default: {"tolwfr": 1.0e-22}.
        dde_tol: a dictionary with a single key defining the type of tolerance used for the DDE calculations and
            its value. Default: {"tolvrs": 1.0e-22}.
        skip_dte_permutations: Since the current version of abinit always performs all the permutations of the
            perturbations, even if only one is asked, if True avoids the creation of inputs that will produce
            duplicated outputs.
        manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
    """
    gs_inp = gs_inp.deepcopy()
    gs_inp.pop_irdvars()

    if ph_tol is None:
        ph_tol = {"tolvrs": 1.0e-22}

    if ddk_tol is None:
        ddk_tol = {"tolwfr": 1.0e-22}

    if dde_tol is None:
        dde_tol = {"tolvrs": 1.0e-22}

    multi = []

    multi_ddk = gs_inp.make_ddk_inputs(tolerance=ddk_tol)
    multi_ddk.add_tags(DDK)
    multi.extend(multi_ddk)
    multi_dde = gs_inp.make_dde_inputs(dde_tol, use_symmetries=False, manager=manager)
    multi_dde.add_tags(DDE)
    multi.extend(multi_dde)

    if use_phonons:
        multi_ph = gs_inp.make_ph_inputs_qpoint([0,0,0], ph_tol, manager=manager)
        multi_ph.add_tags(PH_Q_PERT)
        multi.extend(multi_ph)

    # non-linear calculations do not accept more bands than those in the valence. Set the correct values.
    # Do this as last, so not to interfere with the the generation of the other steps.
    nval = gs_inp.structure.num_valence_electrons(gs_inp.pseudos)
    nval -= gs_inp['charge']
    nband = int(round(nval / 2))
    gs_inp.set_vars(nband=nband)
    gs_inp.pop('nbdbuf', None)
    multi_dte = gs_inp.make_dte_inputs(phonon_pert=use_phonons, skip_permutations=skip_dte_permutations,
                                       manager=manager)
    multi_dte.add_tags(DTE)
    multi.extend(multi_dte)

    multi = MultiDataset.from_inputs(multi)
    multi.add_tags(DFPT)

    return multi


def dfpt_from_gsinput(gs_inp, ph_ngqpt=None, qpoints=None, do_ddk=True, do_dde=True, do_strain=True,
                      do_dte=False, ph_tol=None, ddk_tol=None, dde_tol=None, wfq_tol=None, strain_tol=None,
                      skip_dte_permutations=False, manager=None):
    """
    Returns a list of inputs in the form of a MultiDataset to perform a set of calculations based on DFPT including
    phonons, elastic and non-linear properties. Requires a ground state |AbinitInput| as a starting point.

    It will determine if WFQ files should be calculated for some q points and add the NSCF AbinitInputs to the set.
    The original input is included and the inputs have the following tags, according to their function:
    "scf", "ddk", "dde", "nscf", "ph_q_pert", "strain", "dte", "dfpt".

    N.B. Currently (version 8.8.3) anaddb does not support a DDB containing both 2nd order derivatives with qpoints
    different from gamma AND  3rd oreder derivatives. The calculations could be run, but the global DDB will not
    be directly usable as is.

    Args:
        gs_inp: an |AbinitInput| representing a ground state calculation, likely the SCF performed to get the WFK.
        ph_ngqpt: a list of three integers representing the gamma centered q-point grid used for the calculation.
            If None and qpoint==None the ngkpt value present in the gs_input will be used.
            Incompatible with qpoints.
        qpoints: a list of coordinates of q points in reduced coordinates for which the phonon perturbations will
            be calculated. Incompatible with ph_ngqpt.
        do_ddk: If True, if Gamma is included in the list of qpoints it will add inputs for the calculations of
            the DDK.
        do_dde: If True, if Gamma is included in the list of qpoints it will add inputs for the calculations of
            the DDE. Automatically sets with_ddk=True.
        do_strain: If True inputs for the strain perturbations will be included.
        do_dte: If True inputs for the non-linear perturbations will be included. The phonon non-linear perturbations
            will be included only if a phonon calculation at gamma is present. The caller is responsible for
            adding it.
        ph_tol: a dictionary with a single key defining the type of tolerance used for the phonon calculations and
            its value. Default: {"tolvrs": 1.0e-10}.
        ddk_tol: a dictionary with a single key defining the type of tolerance used for the DDK calculations and
            its value. Default: {"tolwfr": 1.0e-22}.
        dde_tol: a dictionary with a single key defining the type of tolerance used for the DDE calculations and
            its value. Default: {"tolvrs": 1.0e-10}.
        wfq_tol: a dictionary with a single key defining the type of tolerance used for the NSCF calculations of
            the WFQ and its value. Default {"tolwfr": 1.0e-22}.
        strain_tol:  dictionary with a single key defining the type of tolerance used for the strain calculations of
            and its value. Default {"tolvrs": 1.0e-12}.
        skip_dte_permutations: Since the current version of abinit always performs all the permutations of the
            perturbations, even if only one is asked, if True avoids the creation of inputs that will produce
            duplicated outputs.
        manager: |TaskManager| of the task. If None, the manager is initialized from the config file.
    """

    if ph_tol is None:
        ph_tol = {"tolvrs": 1.0e-10}
    if ddk_tol is None:
        ddk_tol = {"tolwfr": 1.0e-22}
    if dde_tol is None:
        dde_tol = {"tolvrs": 1.0e-10}
    if wfq_tol is None:
        wfq_tol = {"tolwfr": 1.0e-22}
    if strain_tol is None:
        strain_tol = {"tolvrs": 1.0e-12}

    if do_dde:
        do_ddk = True

    multi = MultiDataset.from_inputs([gs_inp])
    multi[0].add_tags(SCF)

    do_phonons = ph_ngqpt is not None or qpoints is not None
    has_gamma = False
    if do_phonons:
        multi.extend(phonons_from_gsinput(gs_inp, ph_ngqpt=ph_ngqpt, qpoints=qpoints, with_ddk=False, with_dde=False,
                                          with_bec=False, ph_tol=ph_tol, ddk_tol=ddk_tol, dde_tol=dde_tol,
                                          wfq_tol=wfq_tol, qpoints_to_skip=None, manager=manager))
        has_gamma = ph_ngqpt is not None or any(np.allclose(q, [0, 0, 0]) for q in qpoints)

    if do_ddk:
        multi_ddk = gs_inp.make_ddk_inputs(tolerance=ddk_tol)
        multi_ddk.add_tags(DDK)
        multi.extend(multi_ddk)
    if do_dde:
        multi_dde = gs_inp.make_dde_inputs(dde_tol, use_symmetries=not do_dte, manager=manager)
        multi_dde.add_tags(DDE)
        multi.extend(multi_dde)

    if do_strain:
        multi_strain = gs_inp.make_strain_perts_inputs(tolerance=strain_tol, manager=manager, phonon_pert=False,
                                                       kptopt=2)
        multi_strain.add_tags([DFPT, STRAIN])
        multi.extend(multi_strain)

    if do_dte:
        # non-linear calculations do not accept more bands than those in the valence. Set the correct values.
        nval = gs_inp.structure.num_valence_electrons(gs_inp.pseudos)
        nval -= gs_inp['charge']
        nband = int(round(nval / 2))
        gs_inp_copy = gs_inp.deepcopy()
        gs_inp_copy.set_vars(nband=nband)
        gs_inp_copy.pop('nbdbuf', None)
        multi_dte = gs_inp_copy.make_dte_inputs(phonon_pert=do_phonons and has_gamma,
                                           skip_permutations=skip_dte_permutations, manager=manager)
        multi_dte.add_tags([DTE, DFPT])
        multi.extend(multi_dte)

    return multi

#FIXME if the pseudos are passed as a PseudoTable the whole table will be serialized,
# it would be better to filter on the structure elements
class InputFactory(MSONable):
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
                raise ValueError('An input is required for factory function {0}.'.format(self.factory_function.__name__))
            abiinput = self.factory_function(previous_input, *self.args, **kwargs)
        else:
            abiinput = self.factory_function(*self.args, **kwargs)

        for d in decorators:
            abiinput = d(abiinput)
        abiinput.set_vars(extra_abivars)

        return abiinput

    @pmg_serialize
    def as_dict(self):
        # sanitize to avoid numpy arrays and serialize MSONable objects
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


class HybridScfFromGsFactory(InputFactory):
    factory_function = staticmethod(hybrid_scf_input)


class ScfFactory(InputFactory):
    factory_function = staticmethod(scf_input)
    input_required = False


class ScfForPhononsFactory(InputFactory):
    factory_function = staticmethod(scf_for_phonons)
    input_required = False


class PhononsFromGsFactory(InputFactory):
    factory_function = staticmethod(phonons_from_gsinput)


class PiezoElasticFactory(InputFactory):
    factory_function = staticmethod(scf_piezo_elastic_inputs)
    input_required = False


class PiezoElasticFromGsFactory(InputFactory):
    factory_function = staticmethod(piezo_elastic_inputs_from_gsinput)
