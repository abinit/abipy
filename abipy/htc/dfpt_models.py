"""Pydantic modesl for storing the results of DFPT calculations"""
from __future__ import annotations

from pydantic import Field
from .base_models import AbipyModel
from abipy.dfpt.ddb import DdbFile


class PhononData(AbipyModel):

    ddb_string: str = Field(..., description="DDB string")

    #min_phfreq_ev: float = Field(..., description="DDB string")
    #max_phfreq_ev: float = Field(..., description="DDB string")

    #phbands: PhononBands = Field(..., description="DDB string")
    #phdos: PhononDos = Field(..., description="DDB string")

    @classmethod
    def from_ddb(cls, ddb: DdbFile) -> PhononData:
        kwargs = dict(
            ddb_string=ddb.get_string(),
        )

        #with ddb.anaget_phbst_and_phdos_files() as g:
        #    phbst_file, phdos_file = g[0], g[1]
        #    phbands = phbst_file.phbands
        #    phdos = phdos_file.phdos
        #    kwargs.update(dict(
        #        phbands=phbands,
        #        phdos=phdos,
        #
        #    ))

        return cls(**kwargs)