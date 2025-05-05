"""Tests for ml_phonopy module"""
import os
import numpy as np
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.core.structure import Structure
from abipy.ml.ml_phonopy import MlPhonopyWithDDB, MlPhonopy, MlVZSISAQHAPhonopy


class AbimlTest(AbipyTest):

    def test_mlphonopy_with_ddb(self):
        """Testing MLPhonopyWithDDB"""
        ddb_filepath = os.path.join(abidata.dirpath, "refs", "al_eph", "out_444q_DDB")
        ml = MlPhonopyWithDDB(
                 ddb_filepath,
                 distance=0.01,
                 asr=1,
                 dipdip=0,
                 line_density=10,
                 qppa=None,
                 relax_mode="cell",
                 fmax=0.01,
                 pressure=0,
                 steps=100,
                 optimizer="BFGS",
                 nn_names=["emt"],
                 verbose=1,
                 workdir=None,
                 prefix=None,
                 supercell=None
                )

        assert ml.to_string(verbose=1)
        print(ml)
        ml.run()
        ml.rmtree()

    def test_mlphonopy(self):
        """Testing MlPhonopy."""
        structure = Structure.as_structure(abidata.cif_file("al.cif"))

        ml = MlPhonopy(structure,
                       supercell=[2,2,2],
                       distance=0.01,
                       line_density=10,
                       qppa=None,
                       relax_mode="cell",
                       fmax=0.01,
                       pressure=0.0,
                       steps=100,
                       optimizer="FIRE",
                       nn_names=["emt"],
                       verbose=1,
                       workdir=None,
                       prefix=None
                       )

        assert ml.to_string(verbose=1)
        print(ml)
        ml.run()
        ml.rmtree()

    def test_ml_vzsisa_phonopy(self):
        """Testing MlVZSISAQHAPhonopy."""
        structure = Structure.as_structure(abidata.cif_file("al.cif"))

        bo_vol_scales = [0.96, 0.98, 1, 1.02, 1.04]    # EinfVib4(S)
        ph_vol_scales = [1, 1.02, 1.04]                # EinfVib2(D)
        ph_vol_scales = bo_vol_scales

        ml = MlVZSISAQHAPhonopy(structure,
                               bo_vol_scales,
                               supercell=(2, 2, 2),
                               distance=0.01,
                               line_density=10,
                               qppa=None,
                               relax_mode="cell",
                               fmax=0.01,
                               pressure=0,
                               steps=10,
                               optimizer="BFGS",
                               nn_name="emt",
                               verbose=1,
                               workdir=None,
                               prefix=None
                               )

        assert ml.to_string(verbose=1)
        print(ml)
        ml.run()
        ml.rmtree()
