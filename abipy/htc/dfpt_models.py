"""
Pydantic models storing the results of DFPT calculations
"""
from __future__ import annotations

from pydantic import Field
from abipy.dfpt.ddb import DdbFile
from .base_models import AbipyModel #, MongoConnector, GfsFileDesc, GfsDesc


class PhononData(AbipyModel):

    #min_phfreq_ev: float = Field(..., description="Minimum phonon frequency in eV. Includes 'negative' values")
    #max_phfreq_ev: float = Field(..., description="Maximum phonon frequency in eV.")

    #ddb_gfsd: GfsFileDesc = Field(None, description="Metadata needed to retrieve the DDB file from GridFS.")

    #phbands: PhononBands = Field(..., description="DDB string")
    #phdos: PhononDos = Field(..., description="DDB string")

    #phbands_gfsd: GfsDesc = Field(None, description="Pointer to the PhononBands object in GridFS.")
    #phdos_gfsd: GfsDesc = Field(None, description="Pointer to the PhononDos object in GridFS.")

    @classmethod
    def from_ddb(cls, ddb: DdbFile, mng_connector: MongoConnector) -> PhononData:
        """
        Build the model from a DdbFile `ddb`.
        """
        data = dict(
            #ddb_gfsd=mng_connector.gfs_put_filepath(ddb.filepath)
        )

        # TODO: Add Msonable protocol to phbands and phdos.
        # Should we store PHBST, PHDOS and anaddb.nc or just MSonable objects ?

        #with ddb.anaget_phbst_and_phdos_files() as g:
        #    phbst_file, phdos_file = g[0], g[1]
        #    phbands = phbst_file.phbands
        #    phdos = phdos_file.phdos
        #    data.update(dict(
        #        phbands=phbands,
        #        phdos=phdos,
        #
        #    ))
        #    self.phbands_gfsd = mng_connector.gfs_put_mson(phbst_file.phbands)
        #    self.phdos_gfsd = mng_connector.gfs_put_mson(phdos_file.phdos)

        return cls(**data)
