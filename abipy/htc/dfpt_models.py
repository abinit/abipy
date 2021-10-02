"""
Pydantic models storing the results of DFPT calculations
"""
from __future__ import annotations

from pydantic import Field
from abipy.dfpt.ddb import DdbFile
from .base_models import AbipyModel, MongoConnector, GridFsDesc


class PhononData(AbipyModel):

    #min_phfreq_ev: float = Field(..., description="Minimum phonon frequency in eV. Includes 'negative' values")
    #max_phfreq_ev: float = Field(..., description="Maximum phonon frequency in eV.")

    #phbands: PhononBands = Field(..., description="DDB string")
    #phdos: PhononDos = Field(..., description="DDB string")

    #ddb_string: str = Field(..., description="String with the content of the DDB file")
    ddb_gfsd: GridFsDesc = Field(None, description="Metadata needed to retrieve the DDB file from GridFS.")

    @classmethod
    def from_ddb(cls, ddb: DdbFile, mongo_connector: MongoConnector) -> PhononData:
        """
        Build the model from a DdbFile.
        """
        kwargs = dict(
            #ddb_string=ddb.get_string(),
        )

        kwargs["ddb_gfsd"] = mongo_connector.gridfs_put_filepath(ddb.filepath)

        # TODO: Add Msonable protocol to phbands and phdos.
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
