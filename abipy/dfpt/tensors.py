# coding: utf-8
"""
This modules providess tensors objects extracted from dfpt calculations.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np
from pymatgen.analysis.elasticity.tensors import Tensor, SquareTensor
from pymatgen.core.units import bohr_to_angstrom, eV_to_Ha
from abipy.iotools import ETSF_Reader
from abipy.core.abinit_units import Ha_cmm1


class NLOpticalSusceptibilityTensor(Tensor):
    """
    Subclass of pymatgen.analysis.elasticity.tensors.Tensor containing the
    non-linear optical susceptibility tensor.
    """

    @classmethod
    def from_file(cls, filepath):
        """
        Creates the tensor from a anaddb.nc netcdf file that contains the dchide.
        This requires to run anaddb with tnlflag > 0
        """
        reader = ETSF_Reader(filepath)

        try:
            return cls(reader.read_value("dchide"))
        except Exception as exc:
            import traceback
            msg = traceback.format_exc()
            msg += ("Error while trying to read from file.\n"
                    "Verify that nlflag > 0 in anaddb\n")
            raise ValueError(msg)


class DielectricTensor(object):
    """
    Object describing the frequency dependent dielectric tensor as obtained
    from DFPT calculations. The values are calculated on the fly
    based on the phonon frequencies at gamma and oscillator strengths.
    The first three frequencies would be considered as acoustic modes and
    ignored in the calculation of the quantities. No checks would be performed.
    """

    def __init__(self, phfreqs, oscillator_strength, emacro, structure):
        """
        Args:
             phfreqs: a numpy array containing the 3*(num atoms) phonon frequencies at gamma in atomic units
             oscillator_strength: a complex numpy array with shape (number of phonon modes, 3, 3) in atomic units
             emacro: a numpy array containing the dielectric tensor without frequency dependence
                (at infinite frequency) in atomic units
             structure: a pymatgen Structure of the system considered
        """
        self.phfreqs = phfreqs
        self.oscillator_strength = oscillator_strength
        self.emacro = emacro
        self.structure = structure

    @classmethod
    def from_files(cls, phbst_filepath, anaddbnc_filepath):
        """
        Generates the object from the files that contain the phonon frequencies, oscillator strength and
        static dielectric tensor, i.e. the PHBST and anaddb netcdf files, respectively.
        """
        reader_phbst = ETSF_Reader(phbst_filepath)
        reader_anaddbnc = ETSF_Reader(anaddbnc_filepath)

        qpts = reader_phbst.read_value("qpoints")
        full_phfreqs = reader_phbst.read_value("phfreqs")

        for i, q in enumerate(qpts):
            if np.array_equal(q, [0, 0, 0]):
                break
        else:
            raise ValueError('The PHBST does not containg the frequencies at gamma')

        phfreqs = full_phfreqs[i] * eV_to_Ha

        emacro = reader_anaddbnc.read_value("emacro_cart")

        try:
            oscillator_strength = reader_anaddbnc.read_value("oscillator_strength", cmode="c")
        except Exception as exc:
            import traceback
            msg = traceback.format_exc()
            msg += ("Error while trying to read from file.\n"
                    "Verify that dieflag == 1, 3 or 4 in anaddb\n")
            raise ValueError(msg)

        structure = reader_anaddbnc.read_structure()

        return cls(phfreqs, oscillator_strength, emacro, structure)


    @classmethod
    def from_objects(cls, phbands, anaddbnc):
        """
        Generates the object from the objects PhononBands and AnaddbNcFile
        """
        gamma_index = phbands.qindex([0, 0, 0])

        phfreqs = phbands.phfreqs[gamma_index] * eV_to_Ha

        emacro = anaddbnc.emacro.cartesian_tensor
        oscillator_strength = anaddbnc.oscillator_strength

        return cls(phfreqs, oscillator_strength, emacro, anaddbnc.structure)

    def tensor_at_frequency(self, w, units='Ha'):
        """
        Returns a pymatgen.analysis.elasticity.tensors.SquareTensor object representing
        the dielectric tensor in atomic units at the specified frequency w.
        Args:
            w: frequency
            units: string specifying the units used for the frequency. Accepted values are Ha (default), eV, cm-1
        """

        if units == 'Ha':
            pass
        elif units == 'eV':
            w = w * eV_to_Ha
        elif units == 'cm-1':
            w = w / Ha_cmm1
        else:
            raise ValueError('Unknown units {}'.format(units))

        t = np.zeros((3,3))
        for i in range(3, len(self.phfreqs)):
            t += self.oscillator_strength[i].real/(self.phfreqs[i]**2 - w**2)

        vol = self.structure.volume / bohr_to_angstrom ** 3
        t = 4*np.pi*t/vol

        t += self.emacro

        return SquareTensor(t)