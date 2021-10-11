"""
"""
from __future__ import annotations

from ruamel import yaml

#from abc import ABC, abstractmethod
from typing import Dict, Any, Literal
from pydantic import Field, root_validator
from abipy.abio.abivars import is_abivar, is_anaddb_var
from .base_models import AbipyModel
from .pseudos_models import PseudoSpecs


class Entry(AbipyModel):

    meta_params: Dict[str, Any] = Field(default_factory=dict, description="Dictionary with meta-parameters")

    abivars: Dict[str, Any] = Field(default_factory=dict, description="Dictionary with Abinit variables")

    @root_validator
    def check_abivars(cls, values):
        abivars = values.get("abivars")
        for k in abivars:
            if not is_abivar(k):
                raise ValueError(f"Uknown Variable: `{k}`")
        return values


class Protocol(AbipyModel):
    """
    A protocol contains metaparameters used to generate input files
    """
    #name: str

    description: str = Field(..., description="Human readable description of the protocol.")

    cutoff_stringency: Literal["low", "normal", "high"] = Field(..., description="")

    pseudo_specs: PseudoSpecs = Field(..., description="GS SCF calculation")

    gs_scf: Entry = Field(None, description="GS SCF calculation")

    gs_nscf_kpath: Entry = Field(None, description="NSCF band structure calculations along a k-path")

    gs_nscf_kdos: Entry = Field(None, description="NSCF calculations with k-mesh to compute DOS")

    relax: Entry = Field(None, description="Structural relaxation")

    # validators
    #_normalize_name = validator('name', allow_reuse=True)(normalize)
    #_normalize_name = reuse_validator('name')(normalize)

    #@root_validator
    #def check_pseudos(cls, values):
    #    extra_abivars = values.get("extra_abivars")
    #    return values


yaml_example = """

# Variables common to all the protocols
globals:
  nsppol: 1
  nspden: 1
  nspinor: 1
  #_spin_mode: "unpolarized"
  # default to metallic treatment
  occopt: 7   # Gaussian smearing
  fband: 2.00 # increase the number of bands to > fband * natoms
  nstep: 100
  paral_kgb: 1
  rmm_diis: 1
  expert_user: 1


fast:
    description: "Protocol to relax a structure with low precision at minimal computational cost for testing purposes."

    pseudos:
      repo_name: "ONCVPSP-PBEsol-SR-PDv0.4"
      table_name: "standard"
      cutoff_stringency: "low"

    gs_scf:
      #kpoints_distance: 0.25
      _kppa: 500
      nshiftk: 1
      shiftk: [[0.0, 0.0, 0.0]]
      tolvrs: 1.0e-7
      tsmear: 0.008 # Ha

    gs_nscf_kpath:
      ndivsm: 20
      tolwfr: 1.0e-14

    gs_nscf_kdos:
      #kpoints_distance: 0.25
      _kppa: 1000
      nshiftk: 1
      shiftk: [[0.0, 0.0, 0.0]]
      tolwfr: 1.0e-14

    relax:
      _extends: "gs_scf"
      # default to relaxing atoms only
      optcell: 0    # do not optimize the cell, Abinit default
      ionmov: 22    # optimize ionic positions with L-BFGS
      dilatmx: 1.00 # don't book additional memory for basis set enlargement, Abinit default
      ecutsm: 0.00  # don't smear the energy cutoff, Abinit default

    dfpt:
      _extends: "gs_scf"
      _qppa: 1000
      #ddk_tolwfr:
      #with_becs: yes
      #with_quad: no


moderate:
    description: "Protocol to relax a structure with normal precision at moderate computational cost."
    _extends: "fast"

    pseudos:
      cutoff_stringency: "normal"

    gs_scf:
      #kpoints_distance: 0.20
      _kppa: 1000
      nshiftk: 1
      shiftk: [[0.0, 0.0, 0.0]]
      # tolerances
      tolvrs: 1.0e-9
      tsmear: 0.008  # Ha

    gs_nscf_kpath:
      ndivsm: 20
      tolwfr: 1.0e-18

    gs_nscf_kdos:
      #kpoints_distance: 0.25
      _kppa: 2000
      nshiftk: 1
      shiftk: [[0.0, 0.0, 0.0]]
      tolwfr: 1.0e-18

    relax:
      _extends: "gs_scf"
      # default to relaxing atoms only
      optcell: 0    # do not optimize the cell, Abinit default
      ionmov: 22    # optimize ionic positions with L-BFGS
      dilatmx: 1.00 # don't book additional memory for basis set enlargement, Abinit default
      ecutsm: 0.00  # don't smear the energy cutoff, Abinit default
"""

_registered_tasks = [
        "gs_scf",
        "gs_nscf_kpath",
        "gs_nscf_kdos",
        "relax",
        "dfpt",
]


class ProtocolParser:

    @classmethod
    def from_file(cls, filepath: str) -> ProtocolParser:
        with open(filepath, "rt") as fh:
            document = yaml.safe_load(fh)
            return cls(document)

    @classmethod
    def from_yaml_string(cls, yaml_string: str) -> ProtocolParser:
        document = yaml.safe_load(yaml_string)
        return cls(document)

    def __init__(self, document: dict) -> None:

        # Get global abinit variables first.
        global_abivars = document.pop("global_abivars", {})
        for k in global_abivars:
            if not is_abivar(k):
                raise ValueError(f"Uknown variable: `{k}` found in global_abivars section.")

        # Get pseudo specifications and build model.
        d = document.pop("pseudo_specs", None)
        if d is None:
            raise ValueError(f"Cannot find `pseudo_specs` section in document")
        pseudo_specs = PseudoSpecs.from_repo_table_name(d["repo_name"], d["table_name"])

        protocols = {}
        for k in list(document.keys()):
            d = document.pop(k)
            if not isinstance(d, dict):
                raise TypeError(f"For key {k}: expecting dictionary, got {type(d)}")
            protocols[k] = d

        # Add global variables
        # This means that one can always override the vaue per entry

        for k, d in protocols.items():

           extend_map = {}

           for t in _registered_tasks:
               if t not in d: continue
               if "meta_params" not in d[t]: d[t]["meta_params"] = {}
               if "abivars" not in d[t]: d[t]["abivars"] = {}
               new_dict = global_abivars.copy()
               new_dict.update(d[t]["abivars"])
               d[t]["abivars"] = new_dict
               if "extends" in d[t]:
                   extend_map[t] = d[t]["extends"]

           for t, super_name in extend_map.items():
               super_doc = d[super_name]
               new_metaparams = super_doc["meta_params"].copy()
               new_metaparams.update(d[t]["meta_params"])
               d[t]["meta_params"] = new_metaparams
               new_abivars = super_doc["abivars"].copy()
               new_abivars.update(d[t]["abivars"])
               d[t]["abivars"] = new_abivars

        self.protocols = {k: Protocol(pseudo_specs=pseudo_specs, **data) for k, data in protocols.items()}


    #def get_protocol(self, name):
