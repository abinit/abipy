"""
"""
from __future__ import annotations

import os

from pprint import pformat
from typing import Dict, Any, Literal, List
from ruamel import yaml
from pydantic import Field, root_validator, PrivateAttr
from abipy.core.structure import Structure
from abipy.abio.inputs import AbinitInput
from abipy.abio.abivars import is_abivar
from abipy.abio.factories import gs_input, ebands_input
from abipy.htc.base_models import AbipyModel
from abipy.htc.pseudos_models import PseudoSpecs


class MetaParam(AbipyModel):

    name: str
    exclude_abivars: List[str] = []
    require_abivars: List[str] = []
    exclude_meta: List[str] = []
    require_meta: List[str] = []


_ALL_METAPARAMS = {
    "kppa": MetaParam(name="kppa", exclude_abivars=["nkpt", "ngkpt", "kptrlatt"]),
    "qppa": MetaParam(name="qppa"),
}


class _Specs(AbipyModel):

    meta_params: Dict[str, Any] = Field(
            default_factory=dict,
            description="Dictionary with meta-parameters")

    abivars: Dict[str, Any] = Field(
            default_factory=dict,
            description="Dictionary with extra Abinit variables that will be added to the AbinitInput")

    spin_mode: str = Field(..., description="Spin polarization mode")

    # Private attributes
    #_smearing: str = PrivateAttr()

    def __init__(self, **data):

        abivars = data["abivars"]
        nsppol = abivars.pop("nsppol", 1)
        nspinor = abivars.pop("nspinor", 1)
        nspden = abivars.pop("nspden", 1)

        from pymatgen.io.abinit.abiobjects import _mode2spinvars  #  SpinMode
        for spin_mode, var in _mode2spinvars.items():
            if (nsppol == var.nsppol and nspinor == var.nspinor and nspden == var.nspden):
                break
        else:
            raise ValueError(f"Cannot find spin_mode associated to\n"
                             f"nsppol: {nsppol}, nspinor: {nspinor}, nspden: {nspden}")

        #from pymatgen.io.abinit.abiobjects import Smearing
        #occopt = int(abivars.pop("occopt", 1))
        #tsmear = abivars.pop("tsmear", 0.01)
        #tsmear = float(tsmear)
        #smearing = Smearing(occopt, tsmear)

        super().__init__(spin_mode=spin_mode, **data)
        self.validate_instance()


    def __str__(self) -> str:
        lines = []
        app = lines.append
        app(f"Spec: {self.__class__.__name__}")
        app(f"Doc: {self.__class__.__doc__}")
        app("")
        app(f"meta_params:\n {pformat(self.meta_params)}\n")
        app(f"abivars:\n {pformat(self.abivars)}\n")

        return "\n".join(lines)

    def validate_instance(self):
        cls_name = self.__class__.__name__

        for meta_name in self.meta_params:
            meta = _ALL_METAPARAMS.get(meta_name, None)
            if meta is None:
                raise ValueError(f"In {cls_name}: unknown meta_parameter `{meta_name}`")

            if meta.exclude_abivars:
                for aname in self.abivars:
                    if aname in meta.exclude_abivars:
                        raise ValueError(f"In {cls_name}: metavariable `{meta_name}` "
                                         f"is not compatible with Abinit variable: `{aname}`")

            for aname in meta.require_abivars:
                if aname not in self.abivars:
                    raise ValueError(f"In {cls_name}: metavariable `{meta_name}` "
                                     f"requires Abinit variable: `{aname}`")

            for mp in self.meta_params:
                if mp in meta.exclude_meta:
                    raise ValueError(f"In {cls_name}: metavariable `{meta_name}` "
                                     f"is not compatible with other meta parameter: `{mp}`")

            for mp in meta.require_meta:
                raise ValueError(f"In {cls_name}: metavariable `{meta_name}` requires `{mp}`")

        if "fband" and "nband" in self.abivars:
            raise ValueError("fband and nband cannot be used together.")


    @root_validator
    def check_abivars(cls, values):
        """
        Make sure the names specified in the abivars section are valid.
        """
        for k in values.get("abivars"):
            if not is_abivar(k):
                raise ValueError(f"Unknown Abinit variable: `{k}`")
        return values

    @root_validator
    def check_meta_params(cls, values):
        """
        Check that the meta_params are:

        - registered
        - do not ext

        """
        abivars = values.get("abivars")
        meta_params = values.get("meta_params")

        # Check whether meta_param is supported by this class
        for meta_name in meta_params:
            if meta_name not in cls._supported_meta_params:
                raise ValueError(f"Specs `{cls.__name__}` does not support meta_parameter: `{meta_name}`")

        for meta_name in meta_params:
            meta = _ALL_METAPARAMS.get(meta_name, None)
            if meta is None:
                raise ValueError(f"Unknown meta_parameter: `{meta_name}`")

            #if meta.exclude_abivars and any(vname in meta.exclude_abivars
            #meta.require_abivars
            #meta.exclude_meta
            #meta.require_meta
            #raise ValueError(f"Unknown Abinit variable: `{k}`")

        return values


class GsScfSpecs(_Specs):
    """
    Specifications for ground-state SCF calculations
    """

    _supported_meta_params = [
        "kppa",
    ]



class GsNscfKpathSpecs(_Specs):
    """
    Specifications for ground-state NSCF calculations with k-path (used to compute band structures)
    """

    _supported_meta_params = [
    ]


class GsNscfKdosSpecs(_Specs):
    """
    Specifications for ground-state NSCF calculations with k-mesh (used to compute e-DOS)
    """

    _supported_meta_params = [
        "kppa",
    ]


class RelaxSpecs(_Specs):
    """
    Specifications for Structural relaxations
    """

    _supported_meta_params = [
        "kppa",
    ]


class PhononSpecs(_Specs):
    """
    Specifications for phonon calculations within DFPT
    """

    _supported_meta_params = [
        "qppa",
    ]


KNOWN_SPECS = [
        "gs_scf_specs",
        "gs_nscf_kpath_specs",
        "gs_nscf_kdos_specs",
        "relax_specs",
        #"phonon_specs",
]


# Absolute path to the directory with yaml files.
PROTO_DIRPATH = os.path.join(os.path.dirname(__file__), "protocols")


class Protocol(AbipyModel):
    """
    A protocol contains meta parameters used to generate input files
    """
    info: str = Field(..., description="Human readable string with the description of the protocol.")

    #name: str =  Field(..., description=".")

    #md5: str = Field(..., description="")

    accuracy: Literal["low", "normal", "high"] = Field(..., description="")

    pseudos_specs: PseudoSpecs = Field(..., description="Pseudopotential Specs")

    gs_scf_specs: GsScfSpecs = Field(..., description="Specs for GS SCF calculation")

    gs_nscf_kpath_specs: GsNscfKpathSpecs = Field(
            None,
            description="Specs for NSCF band structure calculations along a k-path")

    gs_nscf_kdos_specs: GsNscfKdosSpecs = Field(
            None,
            description="Specs for NSCF calculations with k-mesh to compute DOS")

    relax_specs: RelaxSpecs = Field(
            None,
            description="Specs for structural relaxations")

    # validators
    #_normalize_name = validator('name', allow_reuse=True)(normalize)
    #_normalize_name = reuse_validator('name')(normalize)

    #@root_validator
    #def check_pseudo_specs(cls, values):
    #    pseudos_specs = values.get("pseudos_specs")
    #    return values

    @classmethod
    def from_file(cls, filepath: str) -> Protocol:
        """Build a protocol object from a Yaml file."""
        with open(filepath, "rt") as fh:
            data = yaml.safe_load(fh)
            return cls.parse_data(data)

    @classmethod
    def from_yaml_string(cls, yaml_string: str) -> Protocol:
        """Build a protocol object from a string in Yaml format."""
        data = yaml.safe_load(yaml_string)
        return cls.parse_data(data)

    @classmethod
    def from_name(cls, name: str) -> Protocol:
        """Build a protocol object from the name of the Yaml file in the AbiPy library."""
        filepath = os.path.join(PROTO_DIRPATH, name)
        return cls.from_file(filepath)

    @classmethod
    def get_all_abipy_protocols(cls) -> List[Protocol]:
        """
        Return list of all registered protocols.
        """
        yaml_files = [f for f in os.listdir(PROTO_DIRPATH) if f.endswith(".yml")]

        protocols = []
        for y in yaml_files:
            proto = cls.from_file(os.path.join(PROTO_DIRPATH, y))
            protocols.append(proto)

        return protocols

    @classmethod
    def parse_data(cls, data: dict) -> Protocol:
        """
        Parse dictionary and return Protocol.
        """
        data = data.copy()

        #name = data.pop("name", None)
        #if name is None:
        #    raise ValueError(f"Cannot find `name` section in data")

        info = data.pop("info", None)
        if info is None:
            raise ValueError(f"Cannot find `info` key in document")

        accuracy = data.pop("accuracy", None)
        if accuracy is None:
            raise ValueError(f"Cannot find `accuracy` section in document")

        # Get global abinit variables first.
        global_abivars = data.pop("global_abivars", {})

        # Get pseudo specifications and build model.
        d = data.pop("pseudos_specs", None)
        if d is None:
            raise ValueError(f"Cannot find `pseudos_specs` section in data")
        pseudos_specs = PseudoSpecs.from_repo_table_name(d["repo_name"], d["table_name"])

        specs = {}
        for spec_name in KNOWN_SPECS:
            specs[spec_name] = data.pop(spec_name, None)

        if data:
            raise ValueError(f"Found unknown sections with keys: `{list(data.keys())}`")

        # Add global variables
        # This means that one can always override the value per entry
        extend_map = {}
        for spec_name, d in specs.items():
            if d is None: continue
            for n in ["meta_params", "abivars"]:
                if n not in d: d[n] = {}
            new_dict = global_abivars.copy()
            new_dict.update(d["abivars"])
            d["abivars"] = new_dict
            extends = d.pop("extends", None)
            if extends:
                extend_map[spec_name] = str(extends)

        # Implement "extends" syntax.
        for spec_name, super_name in extend_map.items():
            super_d = d[super_name]
            #new_metaps = super_d["meta_params"].copy()
            #new_metaps.update(d[spec_name]["meta_params"])
            #d[spec_name]["meta_params"] = new_metaps
            #new_abivars = super_d["abivars"].copy()
            #new_abivars.update(d[spec_name]["abivars"])
            #d[spec_name]["abivars"] = new_abivars

        return cls(info=info, accuracy=accuracy, pseudos_specs=pseudos_specs, **specs)

    def __init__(self, **data):
        super().__init__(**data)
        self.validate_instance()

    def validate_instance(self):
        cls_name = self.__class__.__name__
        #raise ValueError(f"In {cls_name}: unknown meta_parameter `{meta_name}`")

    #def __repr__(self) -> str:
    #    return "<>"

    def __str__(self) -> str:
        lines = [f"info:\n{self.info}"]
        app = lines.append
        app(f"pseudos: {self.pseudos_specs}")

        for spec_name in KNOWN_SPECS:
            spec = getattr(self, spec_name)
            app(f"spec_name: {spec_name}")
            app(f"spec: {spec}")
            app("")

        return "\n".join(lines)

    ####################
    # Factory functions.
    ####################

    def get_gs_scf_input(self, structure: Structure) -> AbinitInput:
        """
        Build and return an input for a GS SCF calculation for the given structure
        """
        pseudos = self.pseudos_specs.get_pseudos()
        scf_specs = self.gs_scf_specs

        scf_inp = gs_input(structure, pseudos,
                           #kppa=scf_specs.kppa,
                           ecut=None, pawecutdg=None, scf_nband=None,
                           accuracy=self.accuracy,
                           spin_mode=scf_specs.spin_mode,
                           smearing="fermi_dirac:0.1 eV",
                           charge=0.0,
                           )

        # Add extra abivars
        scf_inp.set_vars(**scf_specs.abivars)

        return scf_inp

    def get_ebands_input(self, structure: Structure):
        """
        Build and return list of inputs for a GS SCF + NSCF for the given structure
        """
        pseudos = self.pseudos_specs.get_pseudos()

        scf_specs = self.gs_scf_specs
        nscf_specs = self.gs_nscf_kpath_specs
        if nscf_specs is None:
            raise ValueError("band structure calculation requires the specification of `gs_nscf_kpath_specs`")

        multi = ebands_input(structure, pseudos,
                             kppa=self.kppa,
                             nscf_nband=None,
                             ndivsm=self.ndivsm,
                             ecut=None, pawecutdg=None,
                             scf_nband=None,
                             accuracy=self.accuracy,
                             spin_mode=self.spin_mode,
                             smearing=self.smearing,
                             charge=self.charge,
                             #dos_kppa=self.dos_kppa,
                             )

        #nscf_inp = scf_inp.make_ebands_input(ndivsm=15, tolwfr=1e-20, nscf_nband=None, nb_extra=10)
        scf_inp, nscf_inp = multi.split_datasets()

        # Add extra abivars
        scf_inp.set_vars(**self.gs_scf_specs.abivars)
        nscf_inp.set_vars(**self.gs_nscf_specs.abivars)

        return scf_inp, nscf_inp

    #def get_relax_input(self, structure: Structure) -> AbinitInput:
    #    """
    #    Build and return an input to relax the input structure
    #    """
    #    scf_inp = self.get_gs_scf_input(structure)
    #    relax_specs = self.relax_specs
    #    if relax_specs is None:
    #         raise ValueError("structure relaxations require the specifications of `relax_specs`")
    #    relax_input = None
    #    return relax_input


if __name__ == "__main__":
    proto = Protocol.from_file("protocols/NC-PBE-SR-PDv0.4.yml")
    print(proto)

    from abipy.data.ucells import structure_from_ucell
    structures = [structure_from_ucell(name) for name in ("Si",)] # "Si-shifted")]
    for structure in structures:
        scf_inp = proto.get_gs_scf_input(structure)
        print(scf_inp)
