# coding: utf-8
"""
Objects common to the other eph modules.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
import abipy.core.abinit_units as abu

from collections import OrderedDict
from monty.functools import lazy_property
from abipy.electrons.ebands import ElectronsReader

# Phonon frequency in Ha below which e-ph matrix elements are set to zero.
EPH_WTOL = 1e-6

class BaseEphReader(ElectronsReader):
    """
    Provides methods common to the netcdf files produced by the EPH code.
    See Abinit docs for the meaning of the variables.
    """

    @lazy_property
    def ddb_ngqpt(self):
        """Q-Mesh for DDB file."""
        return self.read_value("ddb_ngqpt")

    @lazy_property
    def ngqpt(self):
        """Effective Q-mesh used in to compute integrals (ph_linewidts, e-ph self-energy)."""
        return self.read_value("ngqpt")

    @lazy_property
    def ph_ngqpt(self):
        """Q-mesh for Phonon DOS, interpolated A2F ..."""
        return self.read_value("ph_ngqpt")

    @lazy_property
    def eph_ngqpt_fine(self):
        """Q-mesh for interpolated DFPT potentials"""
        return self.read_value("eph_ngqpt_fine")

    @lazy_property
    def common_eph_params(self):
        """
        Read basic parameters (scalars) from the netcdf files produced by the EPH code and cache them
        """
        od = OrderedDict([
            ("ddb_nqbz", np.prod(self.ddb_ngqpt)),
            ("eph_nqbz_fine", np.prod(self.eph_ngqpt_fine)),
            ("ph_nqbz", np.prod(self.ph_ngqpt)),
        ])

        for vname in ["eph_intmeth", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie"]:
            value = self.read_value(vname)
            if vname in ("eph_intmeth",):
                value = int(value)
            else:
                value = float(value)
            od[vname] = value

        return od


def glr_frohlich(qpoint, becs_cart, epsinf_cart, phdispl_cart_bohr, phfreqs_ha, structure, 
                 qdamp=None, eph_wtol=EPH_WTOL, tol_qnorm=1e-6):
    """
    Compute the long-range part of the e-ph matrix element with the simplified Frohlich model 
    i.e. we include only G = 0 and the <k+q,b1|e^{i(q+G).r}|b2,k> coefficient is replaced by delta_{b1, b2}

    Args:
        qpoint: |Kpoint| object.
        becs_cart: (natom, 3, 3) arrays with Born effective charges in Cartesian coordinates.
        epsinf_cart: (3, 3) array with macroscopic dielectric tensor in Cartesian coordinates.
        phdispl_cart_bohr: (natom3_nu, natom3) complex array with phonon displacement in Cartesian coordinates (Bohr)
        phfreqs_ha: (3 * natom) array with phonon frequencies in Ha.
        structure: |Structure| object.
        qdamp: Exponential damping.
        eph_wtol: Set g to zero below this phonon frequency.
        tol_qnorm: Tolerance of the norm of the q-point.

    Return:
        (natom3) complex array with glr_nu.
    """
    natom = len(structure)
    natom3 = natom * 3
    # 2pi factor already included in frac_coords
    qcs = qpoint.cart_coords
    q_eps_q = np.dot(qcs, np.matmul(epsinf_cart, qcs))
    phdispl_cart_bohr = np.reshape(phdispl_cart_bohr, (natom3, natom, 3))

    xred = structure.frac_coords
    # Acoustic modes are included --> assume BECS fullfill charge neutrality
    glr_nu = np.empty(natom3, dtype=np.complex)
    for nu in range(natom3):
        if phfreqs_ha[nu] < EPH_WTOL or q_eps_q < tol_qnorm: continue
        num = 0.0j
        for iat in range(natom):
            cdd = phdispl_cart_bohr[nu, iat] * np.exp(-2.0j * np.pi * np.dot(qpoint.frac_coords, xred[iat]))
            num += np.dot(qcs, np.matmul(becs_cart[iat], cdd))
        glr_nu[nu] = num / (q_eps_q * np.sqrt(2.0 * phfreqs_ha[nu]))

    fact = 1
    if (qdamp is not None): fact = np.exp(-qpoint.norm**2/(4*qdamp))
    return fact * glr_nu * 4j * np.pi / (structure.volume*abu.Ang_Bohr**3)
