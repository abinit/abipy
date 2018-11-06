# coding: utf-8
"""DDB File."""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import tempfile
import itertools
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from collections import OrderedDict
from six.moves import map, zip
from monty.string import marquee, list_strings
from monty.collections import AttrDict, dict2namedtuple, tree
from monty.functools import lazy_property
from monty.termcolor import cprint
from monty.dev import deprecated
from pymatgen.core.units import eV_to_Ha, bohr_to_angstrom, ang_to_bohr, Energy
from abipy.flowtk import NetcdfReader, AnaddbTask
from abipy.core.mixins import TextFile, Has_Structure, NotebookWriter
from abipy.core.symmetries import AbinitSpaceGroup
from abipy.core.structure import Structure
from abipy.core.kpoints import KpointList, Kpoint
from abipy.iotools import ETSF_Reader
from abipy.tools.numtools import data_from_cplx_mode
from abipy.abio.inputs import AnaddbInput
from abipy.dfpt.phonons import PhononDosPlotter, PhononBandsPlotter
from abipy.dfpt.ifc import InteratomicForceConstants
from abipy.dfpt.elastic import ElasticData
from abipy.core.abinit_units import phfactor_ev2units, phunit_tag
from abipy.tools.plotting import Marker, add_fig_kwargs, get_ax_fig_plt, set_axlims, get_axarray_fig_plt
from abipy.tools import duck
from abipy.tools.tensors import DielectricTensor, ZstarTensor, Stress
from abipy.abio.robots import Robot

try:
    from functools import lru_cache
except ImportError:  # py2k
    from abipy.tools.functools_lru_cache import lru_cache


class DdbError(Exception):
    """Error class raised by DDB."""


class AnaddbError(DdbError):
    """
    Exceptions raised when we try to execute :class:`AnaddbTask` in the :class:`DdbFile` methods

    An `AnaddbError` has a reference to the task and to the :class:`EventsReport` that contains
    the error messages of the run.
    """
    def __init__(self, *args, **kwargs):
        self.task, self.report = kwargs.pop("task"), kwargs.pop("report")
        super(AnaddbError, self).__init__(*args, **kwargs)

    def __str__(self):
        lines = ["\nworkdir = %s" % self.task.workdir]
        app = lines.append

        if self.report.errors:
            app("Found %d errors" % len(self.report.errors))
            lines += [str(err) for err in self.report.errors]

        return "\n".join(lines)


class DdbFile(TextFile, Has_Structure, NotebookWriter):
    """
    This object provides an interface to the DDB_ file produced by ABINIT
    as well as methods to compute phonon band structures, phonon DOS, thermodinamical properties ...

    About the indices (idir, ipert) used by Abinit (Fortran notation):

    * idir in [1, 2, 3] gives the direction (usually reduced direction, cart for strain)
    * ipert in [1, 2, ..., mpert] where mpert = natom + 6

        * ipert in [1, ..., natom] corresponds to atomic perturbations  (reduced dirs)
        * ipert = natom + 1 gives d/dk  (reduced dirs)
        * ipert = natom + 2 gives the electric field
        * ipert = natom + 3 gives the uniaxial stress (cartesian dirs)
        * ipert = natom + 4 gives the shear strees.   (cartesian dirs)

    .. rubric:: Inheritance
    .. inheritance-diagram:: DdbFile
    """
    Error = DdbError
    AnaddbError = AnaddbError

    @classmethod
    def from_file(cls, filepath):
        """Needed for the :class:`TextFile` abstract interface."""
        return cls(filepath)

    @classmethod
    def from_mpid(cls, material_id, api_key=None, endpoint=None):
        """
        Fetch DDB file corresponding to a materials project ``material_id``,
        save it to temporary file and return new DdbFile object.

        Raises: MPRestError if DDB file is not available

        Args:
            material_id (str): Materials Project material_id (e.g., mp-1234).
            api_key (str): A String API key for accessing the MaterialsProject REST interface.
                If None, the code will check if there is a `PMG_MAPI_KEY` in your .pmgrc.yaml.
            endpoint (str): Url of endpoint to access the MaterialsProject REST interface.
                Defaults to the standard Materials Project REST address
        """
        from abipy.core import restapi
        with restapi.get_mprester(api_key=api_key, endpoint=endpoint) as rest:
            ddb_string = rest._make_request("/materials/%s/abinit_ddb" % material_id)

        _, tmpfile = tempfile.mkstemp(prefix=material_id, suffix='_DDB')
        with open(tmpfile, "wt") as fh:
            fh.write(ddb_string)

        return cls(tmpfile)

    @classmethod
    def as_ddb(cls, obj):
        """
        Return an instance of |DdbFile| from a generic object `obj`.
        Accepts: DdbFile or filepath
        """
        return obj if isinstance(obj, cls) else cls.from_file(obj)

    def __init__(self, filepath):
        super(DdbFile, self).__init__(filepath)

        self._header = self._parse_header()

        self._structure = Structure.from_abivars(**self.header)
        # Add AbinitSpacegroup (needed in guessed_ngkpt)
        # FIXME: kptopt is not reported in the header --> has_timerev is always set to True
        spgid, has_timerev, h = 0, True, self.header
        self._structure.set_abi_spacegroup(AbinitSpaceGroup(spgid, h.symrel, h.tnons, h.symafm, has_timerev))

        # Add forces to structure.
        if self.cart_forces is not None:
            self._structure.add_site_property("cartesian_forces", self.cart_forces)

        frac_coords = self._read_qpoints()
        self._qpoints = KpointList(self.structure.lattice.reciprocal_lattice, frac_coords, weights=None, names=None)

    def __str__(self):
        """String representation."""
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []
        app, extend = lines.append, lines.extend

        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(marquee("DDB Info", mark="="))
        app("")
        app("Number of q-points in DDB: %d" % len(self.qpoints))
        app("guessed_ngqpt: %s (guess for the q-mesh divisions made by AbiPy)" % self.guessed_ngqpt)
        #if verbose:
        h = self.header
        #app("Important parameters extracted from the header:")
        app("ecut = %f, ecutsm = %f, nkpt = %d, nsym = %d, usepaw = %d" % (h.ecut, h.ecutsm, h.nkpt, h.nsym, h.usepaw))
        app("nsppol %d, nspinor %d, nspden %d, ixc = %d, occopt = %d, tsmear = %f" % (
            h.nsppol, h.nspinor, h.nspden, h.ixc, h.occopt, h.tsmear))
        app("")

        app("Has total energy: %s, Has forces: %s" % (
            self.total_energy is not None, self.cart_forces is not None))
        if self.total_energy is not None:
            app("Total energy: %s [eV]" % self.total_energy)
        #app("Has forces: %s" % (
        #if self.cart_forces is not None:
        #    app("Cartesian forces (eV/Ang):\n%s" % (self.cart_forces))
        #    app("")
        if self.cart_stress_tensor is not None:
            app("")
            app("Cartesian stress tensor in GPa with pressure %.3e (GPa):\n%s" % (
                - self.cart_stress_tensor.trace() / 3, self.cart_stress_tensor))
        else:
            app("Has stress tensor: %s" % (self.cart_stress_tensor is not None))
        app("")
        app("Has (at least one) atomic pertubation: %s" % self.has_at_least_one_atomic_perturbation())
        #app("Has (at least one) electric-field perturbation: %s" % self.has_epsinf_terms(select="at_least_one"))
        #app("Has (all) electric-field perturbation: %s" % self.has_epsinf_terms(select="all"))
        app("Has (at least one diagonal) electric-field perturbation: %s" % self.has_epsinf_terms(select="at_least_one_diagoterm"))
        app("Has (at least one) Born effective charge: %s" % self.has_bec_terms(select="at_least_one"))
        app("Has (all) strain terms: %s" % self.has_strain_terms(select="all"))
        app("Has (all) internal strain terms: %s" % self.has_internalstrain_terms(select="all"))
        app("Has (all) piezoelectric terms: %s" % self.has_piezoelectric_terms(select="all"))

        if verbose:
            # Print q-points
            app(self.qpoints.to_string(verbose=verbose, title="Q-points in DDB"))

        if verbose > 1:
            # Print full header.
            from pprint import pformat
            app(marquee("DDB Header", mark="="))
            app(pformat(self.header))

        return "\n".join(lines)

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    @property
    def natom(self):
        """Number of atoms in structure."""
        return len(self.structure)

    @property
    def version(self):
        """DDB Version number (integer)."""
        return self.header["version"]

    @property
    def header(self):
        """
        Dictionary with the values reported in the header section.
        Use e.g. ``ddb.header.ecut`` to access its values
        """
        return self._header

    def _parse_header(self):
        """Parse the header sections. Returns |AttrDict| dictionary."""
        #ixc         7
        #kpt  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00
        #     0.25000000000000D+00  0.00000000000000D+00  0.00000000000000D+00
        self.seek(0)
        keyvals = []
        header_lines = []
        for i, line in enumerate(self):
            header_lines.append(line.rstrip())
            line = line.strip()
            if not line: continue
            if "Version" in line:
                # +DDB, Version number    100401
                version = int(line.split()[-1])

            if line in ("Description of the potentials (KB energies)",
                        "No information on the potentials yet",
                        "Description of the PAW dataset(s)"):
                # Skip section with psps info.
                break

            # header starts here
            if i >= 6:
                # Python does not support exp format with D
                line = line.replace("D+", "E+").replace("D-", "E-")
                tokens = line.split()
                key = None
                try:
                    try:
                        float(tokens[0])
                        parse = float if "." in tokens[0] else int
                        keyvals[-1][1].extend(list(map(parse, tokens)))
                    except ValueError:
                        # We have a new key
                        key = tokens.pop(0)
                        parse = float if "." in tokens[0] else int
                        keyvals.append((key, list(map(parse, tokens))))
                except Exception as exc:
                    raise RuntimeError("Exception:\n%s\nwhile parsing ddb header line:\n%s" %
                                        (str(exc), line))

        # add the potential information
        for line in self:
            if "Database of total energy derivatives" in line:
                break
            header_lines.append(line.rstrip())

        h = AttrDict(version=version, lines=header_lines)
        for key, value in keyvals:
            if len(value) == 1: value = value[0]
            h[key] = value

        # Convert to array. Note that znucl is converted into integer
        # to avoid problems with pymatgen routines that expect integer Z
        # This of course will break any code for alchemical mixing.
        arrays = {
            "acell": dict(shape=(3, ), dtype=np.double),
            "amu": dict(shape=(h.ntypat, ), dtype=np.double),
            "kpt": dict(shape=(h.nkpt, 3), dtype=np.double),
            "ngfft": dict(shape=(3, ), dtype=np.int),
            # This is problematic because not all occupation factors are written
            #"occ": dict(shape=(h.nsppol, h.nkpt, h.nband), dtype=np.double),
            "rprim": dict(shape=(3, 3), dtype=np.double),
            "spinat": dict(shape=(h.natom, 3), dtype=np.double),
            "symrel": dict(shape=(h.nsym, 3, 3), dtype=np.int),
            "tnons": dict(shape=(h.nsym, 3), dtype=np.double),
            "xred":  dict(shape=(h.natom, 3), dtype=np.double),
            # In principle these two quantities are double but here we convert to int
            # Alchemical mixing is therefore ignored.
            "znucl": dict(shape=(h.ntypat,), dtype=np.int),
            "zion": dict(shape=(h.ntypat,), dtype=np.int),
            "symafm": dict(shape=(h.nsym,), dtype=np.int),
            "wtk": dict(shape=(h.nkpt,), dtype=np.double),
        }

        for k, ainfo in arrays.items():
            try:
                h[k] = np.reshape(np.array(h[k], dtype=ainfo["dtype"]), ainfo["shape"])
            except Exception as exc:
                print("While Trying to reshape", k)
                raise exc

        # Transpose symrel because Abinit write matrices by columns.
        h.symrel = np.array([s.T for s in h.symrel])

        return h

    def _read_qpoints(self):
        """Read the list q-points from the DDB file. Returns |numpy-array|."""
        # 2nd derivatives (non-stat.)  - # elements :      36
        # qpt  2.50000000E-01  0.00000000E+00  0.00000000E+00   1.0

        # Since there are multiple occurrences of qpt in the DDB file
        # we use seen to remove duplicates.
        self.seek(0)
        tokens, seen = [], set()

        for line in self:
            line = line.strip()
            if line.startswith("qpt") and line not in seen:
                seen.add(line)
                tokens.append(line.replace("qpt", ""))

        qpoints, weights = [], []
        for tok in tokens:
            nums = list(map(float, tok.split()))
            qpoints.append(nums[:3])
            weights.append(nums[3])

        return np.reshape(qpoints, (-1, 3))

    @lazy_property
    def computed_dynmat(self):
        """
        :class:`OrderedDict` mapping q-point object to --> pandas Dataframe.
        The |pandas-DataFrame| contains the columns: "idir1", "ipert1", "idir2", "ipert2", "cvalue"
        and (idir1, ipert1, idir2, ipert2) as index.

        .. note::

            The indices follow the Abinit (Fortran) notation so they start at 1.
        """
        # TODO: Create mapping [(idir1, ipert1), (idir2, ipert2)] --> element
        df_columns = "idir1 ipert1 idir2 ipert2 cvalue".split()

        dynmat = OrderedDict()
        for block in self.blocks:
            # skip the blocks that are not related to second order derivatives
            first_line = block["data"][0].strip()
            if not first_line.startswith("2nd derivatives"):
                continue

            # Build q-point object.
            qpt = Kpoint(frac_coords=block["qpt"], lattice=self.structure.reciprocal_lattice, weight=None, name=None)

            # Build pandas dataframe with df_columns and (idir1, ipert1, idir2, ipert2) as index.
            # Each line in data represents an element of the dynamical matric
            # idir1 ipert1 idir2 ipert2 re_D im_D
            df_rows, df_index = [], []
            for line in block["data"]:
                line = line.strip()
                if line.startswith("2nd derivatives") or line.startswith("qpt"):
                    continue
                try:
                    toks = line.split()
                    idir1, ipert1 = p1 = (int(toks[0]), int(toks[1]))
                    idir2, ipert2 = p2 = (int(toks[2]), int(toks[3]))
                    toks[4] = toks[4].replace("D", "E")
                    toks[5] = toks[5].replace("D", "E")
                    cvalue = float(toks[4]) + 1j*float(toks[5])
                except Exception as exc:
                    cprint("exception while parsing line: %s" % line, "red")
                    raise exc

                df_index.append(p1 + p2)
                df_rows.append(dict(idir1=idir1, ipert1=ipert1, idir2=idir2, ipert2=ipert2, cvalue=cvalue))

            dynmat[qpt] = pd.DataFrame(df_rows, index=df_index, columns=df_columns)

        return dynmat

    @lazy_property
    def blocks(self):
        """
        DDB blocks. List of dictionaries, Each dictionary contains the following keys.
        "qpt" with the reduced coordinates of the q-point.
        "data" that is a list of strings with the entries of the dynamical matrix for this q-point.
        """
        return self._read_blocks()

    def _read_blocks(self):
        # skip until the beginning of the db
        self.seek(0)
        while "Number of data blocks" not in self._file.readline():
            pass

        blocks = []
        block_lines = []
        qpt = None
        dord = None

        for line in self:
            # skip empty lines
            if line.isspace():
                continue

            if "List of bloks and their characteristics" in line:
                # add last block when we reach the last part of the file.
                # This line is present only if DDB has been produced by mrgddb
                if block_lines:
                    blocks.append({"data": block_lines, "qpt": qpt, "dord": dord})
                    block_lines = []
                    qpt = None
                break

            # Don't use lstring because we may reuse block_lines to write new DDB.
            line = line.rstrip()

            # new block --> detect order
            if "# elements" in line:
                if block_lines:
                    blocks.append({"data": block_lines, "qpt": qpt, "dord": dord})

                tokens = line.split()
                num_elements = int(tokens[-1])
                s = " ".join(tokens[:2])
                dord = {"Total energy": 0,
                        "1st derivatives": 1,
                        "2nd derivatives": 2,
                        "3rd derivatives": 3}.get(s, None)
                if dord is None:
                    raise RuntimeError("Cannot detect derivative order from string: `%s`" % s)

                block_lines = []
                qpt = None

            block_lines.append(line)
            if "qpt" in line:
                qpt = list(map(float, line.split()[1:4]))

        if block_lines:
            blocks.append({"data": block_lines, "qpt": qpt, "dord": dord})

        return blocks

    @property
    def qpoints(self):
        """|KpointList| object with the list of q-points in reduced coordinates."""
        return self._qpoints

    def has_qpoint(self, qpoint):
        """True if the DDB file contains this q-point."""
        #qpoint = Kpoint.as_kpoint(qpoint, self.structure.reciprocal_lattice)
        try:
            self.qpoints.index(qpoint)
            return True
        except ValueError:
            return False

    def qindex(self, qpoint):
        """
        The index of the q-point in the internal list of q-points.
        Accepts: |Kpoint| instance or integer.
        """
        if duck.is_intlike(qpoint):
            return int(qpoint)
        else:
            return self.qpoints.index(qpoint)

    @lazy_property
    def guessed_ngqpt(self):
        """
        Guess for the q-mesh divisions (ngqpt) inferred from the list of
        q-points found in the DDB file.

        .. warning::

            The mesh may not be correct if the DDB file contains points belonging
            to different meshes and/or the Q-mesh is shifted.
        """
        return self._guess_ngqpt()

    def _guess_ngqpt(self):
        """
        This function tries to figure out the value of ngqpt from the list of
        points reported in the DDB file.
        """
        if not self.qpoints: return None
        # Build the union of the stars of the q-points.
        all_qpoints = np.empty((len(self.qpoints) * len(self.structure.abi_spacegroup), 3))
        count = 0
        for qpoint in self.qpoints:
            for op in self.structure.abi_spacegroup:
                all_qpoints[count] = op.rotate_k(qpoint.frac_coords, wrap_tows=False)
                count += 1

        # Replace zeros with np.inf
        for q in all_qpoints:
            q[q == 0] = np.inf

        # Compute the minimum of the fractional coordinates along the 3 directions and invert
        #print(all_qpoints)
        smalls = np.abs(all_qpoints).min(axis=0)
        smalls[smalls == 0] = 1
        ngqpt = np.rint(1 / smalls)
        ngqpt[ngqpt == 0] = 1

        return np.array(ngqpt, dtype=np.int)

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        names = ("nkpt", "nsppol", "ecut", "tsmear", "occopt", "ixc", "nband", "usepaw")
        od = OrderedDict()
        for k in names:
            od[k] = self.header[k]
        return od

    def _add_params(self, obj):
        """Add params (meta variable) to object ``obj``. Usually a phonon bands or phonon dos object."""
        if not hasattr(obj, "params"):
            raise TypeError("object %s does not have `params` attribute" % type(obj))
        obj.params.update(self.params)

    @lazy_property
    def total_energy(self):
        """
        Total energy in eV. None if not available.
        """
        for block in self.blocks:
            if block["dord"] == 0:
                ene_ha = float(block["data"][1].split()[0].replace("D", "E"))
                return Energy(ene_ha, "Ha").to("eV")
        return None

    @lazy_property
    def cart_forces(self):
        """
        Cartesian forces in eV / Ang
        None if not available i.e. if the GS DDB has not been merged.
        """
        for block in self.blocks:
            if block["dord"] != 1: continue
            natom = len(self.structure)
            fred = np.empty((natom, 3))
            for line in block["data"][1:]:
                idir, ipert, fval = line.split()[:3]
                # F --> C
                idir, ipert = int(idir) - 1, int(ipert) - 1
                if ipert < natom:
                    fred[ipert, idir] = float(fval.replace("D", "E"))

            # Fred stores d(etotal)/d(xred)
            # this array has *not* been corrected by enforcing
            # the translational symmetry, namely that the sum of force
            # on all atoms is not necessarly zero.
            # Compute fcart using same code as in fred2fcart.
            # Note conversion to cartesian coordinates (bohr) AND
            # negation to make a force out of a gradient.
            gprimd = self.structure.reciprocal_lattice.matrix / (2 * np.pi) * abu.Bohr_Ang
            #fcart = - np.matmul(fred, gprimd)
            fcart = - np.matmul(fred, gprimd.T)
            # Subtract off average force from each force component
            favg = fcart.sum(axis=0) / len(self.structure)
            fcart -= favg

            return fcart * abu.Ha_eV / abu.Bohr_Ang

        return None

    @lazy_property
    def cart_stress_tensor(self):
        """
        |Stress| tensor in cartesian coordinates (GPa units).
        None if not available.
        """
        for block in self.blocks:
            if block["dord"] != 1: continue
            svoigt = np.empty(6)
            # Abinit stress is in cart coords and Ha/Bohr**3
            # Map (idir, ipert) --> voigt
            uniax, shear = len(self.structure) + 3, len(self.structure) + 4
            dirper2voigt = {
                (1, uniax): 0,
                (2, uniax): 1,
                (3, uniax): 2,
                (1, shear): 3,
                (2, shear): 4,
                (3, shear): 5}

            for line in block["data"][1:]:
                idir, ipert, fval = line.split()[:3]
                idp = int(idir), int(ipert)
                if idp in dirper2voigt:
                    svoigt[dirper2voigt[idp]] = float(fval.replace("D", "E"))

            # Convert from Ha/Bohr^3 to GPa
            return Stress.from_voigt(svoigt * abu.HaBohr3_GPa)

        return None

    def has_lo_to_data(self, select="at_least_one"):
        """
        True if the DDB file contains the data required to compute the LO-TO splitting.
        """
        return self.has_epsinf_terms(select=select) and self.has_bec_terms(select=select)

    @lru_cache(typed=True)
    def has_at_least_one_atomic_perturbation(self, qpt=None):
        """
        True if the DDB file contains info on (at least one) atomic perturbation.
        If the coordinates of a q point are provided only the specified qpt will be considered.
        """
        natom = len(self.structure)
        ap_list = list(itertools.product(range(1, 4), range(1, natom + 1)))

        for qpt_dm, df in self.computed_dynmat.items():
            if qpt is not None and qpt_dm != qpt: continue

            index_set = set(df.index)
            for p1 in ap_list:
                for p2 in ap_list:
                    p12 = p1 + p2
                    if p12 in index_set: return True

        return False

    @lru_cache(typed=True)
    def has_epsinf_terms(self, select="at_least_one"):
        """
        True if the DDB file contains info on the electric-field perturbation.

        Args:
            select: Possible values in ["at_least_one", "at_least_one_diagoterm", "all"]
                If select == "at_least_one", we check if there's at least one entry associated to the electric field.
                and we assume that anaddb will be able to reconstruct the full tensor by symmetry.
		"at_least_one_diagoterm" is similar but it only checks for the presence of one diagonal term.
                If select == "all", all tensor components must be present in the DDB file.
        """
        gamma = Kpoint.gamma(self.structure.reciprocal_lattice)
        if gamma not in self.computed_dynmat:
            return False

        index_set = set(self.computed_dynmat[gamma].index)

        natom = len(self.structure)
        ep_list = list(itertools.product(range(1, 4), [natom + 2]))
        for p1 in ep_list:
            for p2 in ep_list:
                p12 = p1 + p2
                p21 = p2 + p1
                if select == "at_least_one":
                    if p12 in index_set: return True
                elif select == "at_least_one_diagoterm":
                    if p12 == p21 and p12 in index_set:
                        return True
                elif select == "all":
                    if p12 not in index_set and p21 not in index_set:
                        return False
                else:
                    raise ValueError("Wrong select %s" % str(select))

        return False if select in ("at_least_one", "at_least_one_diagoterm") else True

    @deprecated(message="has_emacro_terms is deprecated and will be removed in abipy 0.8, use has_epsinf_terms")
    def has_emacro_terms(self, **kwargs):
        return self.has_epsinf_terms(**kwargs)

    @lru_cache(typed=True)
    def has_bec_terms(self, select="at_least_one"):
        """
        True if the DDB file contains info on the Born effective charges.

        Args:
            select: Possible values in ["at_least_one", "all"]
                By default, we check if there's at least one entry associated to atomic displacement
                and electric field and we assume that anaddb will be able to reconstruct the full tensor by symmetry.
                If select == "all", all bec components must be present in the DDB file.
        """
        gamma = Kpoint.gamma(self.structure.reciprocal_lattice)
        if gamma not in self.computed_dynmat:
            return False
        index_set = set(self.computed_dynmat[gamma].index)
        natom = len(self.structure)
        ep_list = list(itertools.product(range(1, 4), [natom + 2]))
        ap_list = list(itertools.product(range(1, 4), range(1, natom + 1)))

        for ap1 in ap_list:
            for ep2 in ep_list:
                p12 = ap1 + ep2
                p21 = ep2 + ap1
                if select == "at_least_one":
                    if p12 in index_set: return True
                elif select == "all":
                    if p12 not in index_set and p21 not in index_set:
                        return False
                else:
                    raise ValueError("Wrong select %s" % str(select))

        return False if select == "at_least_one" else True

    @lru_cache(typed=True)
    def has_strain_terms(self, select="all"):
        """
        True if the DDB file contains info on the (clamped-ion) strain perturbation
        (i.e. 2nd order derivatives wrt strain)

        Args:
            select: Possible values in ["at_least_one", "all"]
                If select == "at_least_one", we check if there's at least one entry associated to the strain.
                and we assume that anaddb will be able to reconstruct the full tensor by symmetry.
                If select == "all", all tensor components must be present in the DDB file.

        .. note::

            As anaddb is not yet able to reconstruct the strain terms by symmetry,
            the default value for select is "all"
        """
        gamma = Kpoint.gamma(self.structure.reciprocal_lattice)
        if gamma not in self.computed_dynmat:
            return False

        index_set = set(self.computed_dynmat[gamma].index)

        natom = len(self.structure)
        sp_list = list(itertools.product(range(1, 4), [natom + 3, natom + 4]))
        for p1 in sp_list:
            for p2 in sp_list:
                p12 = p1 + p2
                p21 = p2 + p1
                if select == "at_least_one":
                    if p12 in index_set: return True
                elif select == "all":
                    if p12 not in index_set and p21 not in index_set:
                        #print("p12", p12, "not in index_set")
                        return False
                else:
                    raise ValueError("Wrong select %s" % str(select))

        return False if select == "at_least_one" else True

    @lru_cache(typed=True)
    def has_internalstrain_terms(self, select="all"):
        """
        True if the DDB file contains internal strain terms
        i.e "off-diagonal" 2nd order derivatives wrt (strain, atomic displacement)

        Args:
            select: Possible values in ["at_least_one", "all"]
                If select == "at_least_one", we check if there's at least one entry associated to the strain.
                and we assume that anaddb will be able to reconstruct the full tensor by symmetry.
                If select == "all", all tensor components must be present in the DDB file.

        .. note::

            As anaddb is not yet able to reconstruct the strain terms by symmetry,
            the default value for select is "all"
        """
        gamma = Kpoint.gamma(self.structure.reciprocal_lattice)
        if gamma not in self.computed_dynmat:
            return False

        index_set = set(self.computed_dynmat[gamma].index)

        natom = len(self.structure)
        sp_list = list(itertools.product(range(1, 4), [natom + 3, natom + 4]))
        ap_list = list(itertools.product(range(1, 4), range(1, natom + 1)))
        for p1 in sp_list:
            for p2 in ap_list:
                p12 = p1 + p2
                p21 = p2 + p1
                if select == "at_least_one":
                    if p12 in index_set: return True
                elif select == "all":
                    if p12 not in index_set and p21 not in index_set:
                        #print("p12", p12, "non in index")
                        return False
                else:
                    raise ValueError("Wrong select %s" % str(select))

        return False if select == "at_least_one" else True

    @lru_cache(typed=True)
    def has_piezoelectric_terms(self, select="all"):
        """
        True if the DDB file contains piezoelectric terms
        i.e "off-diagonal" 2nd order derivatives wrt (electric_field, strain)

        Args:
            select: Possible values in ["at_least_one", "all"]
                If select == "at_least_one", we check if there's at least one entry associated to the strain.
                and we assume that anaddb will be able to reconstruct the full tensor by symmetry.
                If select == "all", all tensor components must be present in the DDB file.

        .. note::

            As anaddb is not yet able to reconstruct the (strain, electric) terms by symmetry,
            the default value for select is "all"
        """
        gamma = Kpoint.gamma(self.structure.reciprocal_lattice)
        if gamma not in self.computed_dynmat:
            return False

        index_set = set(self.computed_dynmat[gamma].index)

        natom = len(self.structure)
        sp_list = list(itertools.product(range(1, 4), [natom + 3, natom + 4]))
        ep_list = list(itertools.product(range(1, 4), [natom + 2]))
        for p1 in sp_list:
            for p2 in ep_list:
                p12 = p1 + p2
                p21 = p2 + p1
                if select == "at_least_one":
                    if p12 in index_set: return True
                elif select == "all":
                    if p12 not in index_set and p21 not in index_set: return False
                else:
                    raise ValueError("Wrong select %s" % str(select))

        return False if select == "at_least_one" else True

    def view_phononwebsite(self, browser=None, verbose=0, dryrun=False, **kwargs):
        """
        Invoke anaddb to compute phonon bands.
        Produce JSON_ file that can be parsed from the phononwebsite_ and open it in ``browser``.

        Args:
            browser: Open webpage in ``browser``. Use default $BROWSER if None.
            verbose: Verbosity level
            dryrun: Activate dryrun mode for unit testing purposes.
            kwargs: Passed to anaget_phbst_and_phdos_files

        Return: Exit status
        """
        # Call anaddb to get phonon bands.
        if "nqsmall" not in kwargs: kwargs["nqsmall"] = 0
        phbst_file, phdos_file = self.anaget_phbst_and_phdos_files(**kwargs)
        phbands = phbst_file.phbands
        phbst_file.close()

        return phbands.view_phononwebsite(browser=browser, verbose=verbose, dryrun=dryrun)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.structure.plot(show=False)
        yield self.qpoints.plot(show=False)
        yield self.structure.plot_bz(show=False)

    def anaget_phmodes_at_qpoint(self, qpoint=None, asr=2, chneut=1, dipdip=1, workdir=None, mpi_procs=1,
                                 manager=None, verbose=0, lo_to_splitting=False, spell_check=True,
                                 directions=None, anaddb_kwargs=None):
        """
        Execute anaddb to compute phonon modes at the given q-point (without LO-TO splitting)

        Args:
            qpoint: Reduced coordinates of the qpoint where phonon modes are computed.
            asr, chneut, dipdip: Anaddb input variable. See official documentation.
            workdir: Working directory. If None, a temporary directory is created.
            mpi_procs: Number of MPI processes to use.
            manager: |TaskManager| object. If None, the object is initialized from the configuration file
            verbose: verbosity level. Set it to a value > 0 to get more information.
            lo_to_splitting: Allowed values are [True, False, "automatic"]. Defaults to False
                If True the LO-TO splitting will be calculated if qpoint == Gamma and the non_anal_directions
                non_anal_phfreqs attributes will be addeded to the phonon band structure.
                "automatic" activates LO-TO if the DDB file contains the dielectric tensor and Born effective charges.
            directions: list of 3D directions along which the LO-TO splitting will be calculated. If None the three
                cartesian direction will be used.
            anaddb_kwargs: additional kwargs for anaddb.

        Return: |PhononBands| object.
        """
        if qpoint is None:
            qpoint = self.qpoints[0]
            if len(self.qpoints) != 1:
                raise ValueError("%s contains %s qpoints and the choice is ambiguous.\n"
                                 "Please specify the qpoint." % (self, len(self.qpoints)))

        # Check if qpoint is in the DDB.
        try:
            iq = self.qindex(qpoint)
        except:
            raise ValueError("input qpoint %s not in %s.\nddb.qpoints:\n%s" % (
                qpoint, self.filepath, self.qpoints))

        qpoint = self.qpoints[iq]

        if lo_to_splitting == "automatic":
            lo_to_splitting = self.has_lo_to_data() and qpoint.is_gamma() and dipdip != 0

        if lo_to_splitting and qpoint.is_gamma() and not self.has_lo_to_data():
            cprint("lo_to_splitting set to True but Eps_inf and Becs are not available in DDB %s:" % self.filepath)

        inp = AnaddbInput.modes_at_qpoint(self.structure, qpoint, asr=asr, chneut=chneut, dipdip=dipdip,
                                          lo_to_splitting=lo_to_splitting, directions=directions,
                                          anaddb_kwargs=anaddb_kwargs, spell_check=spell_check)

        task = self._run_anaddb_task(inp, mpi_procs, workdir, manager, verbose)

        with task.open_phbst() as ncfile:
            if lo_to_splitting and qpoint.is_gamma():
                ncfile.phbands.read_non_anal_from_file(os.path.join(task.workdir, "anaddb.nc"))

            return ncfile.phbands

    def anaget_phbst_and_phdos_files(self, nqsmall=10, qppa=None, ndivsm=20, line_density=None, asr=2, chneut=1, dipdip=1,
                                     dos_method="tetra", lo_to_splitting="automatic", ngqpt=None, qptbounds=None,
                                     anaddb_kwargs=None, verbose=0, spell_check=True,
                                     mpi_procs=1, workdir=None, manager=None):
        """
        Execute anaddb to compute the phonon band structure and the phonon DOS

        Args:
            nqsmall: Defines the homogeneous q-mesh used for the DOS. Gives the number of divisions
                used to sample the smallest lattice vector. If 0, DOS is not computed and
                (phbst, None) is returned.
            qppa: Defines the homogeneous q-mesh used for the DOS in units of q-points per reciproval atom.
                Overrides nqsmall.
            ndivsm: Number of division used for the smallest segment of the q-path.
            line_density: Defines the a density of k-points per reciprocal atom to plot the phonon dispersion.
                Overrides ndivsm.
            asr, chneut, dipdip: Anaddb input variable. See official documentation.
            dos_method: Technique for DOS computation in  Possible choices: "tetra", "gaussian" or "gaussian:0.001 eV".
                In the later case, the value 0.001 eV is used as gaussian broadening.
            lo_to_splitting: Allowed values are [True, False, "automatic"]. Defaults to "automatic"
                If True the LO-TO splitting will be calculated and the non_anal_directions
                and the non_anal_phfreqs attributes will be addeded to the phonon band structure.
                "automatic" activates LO-TO if the DDB file contains the dielectric tensor and Born effective charges.
            ngqpt: Number of divisions for the q-mesh in the DDB file. Auto-detected if None (default).
            qptbounds: Boundaries of the path. If None, the path is generated from an internal database
                depending on the input structure.
            anaddb_kwargs: additional kwargs for anaddb.
            verbose: verbosity level. Set it to a value > 0 to get more information.
            mpi_procs: Number of MPI processes to use.
            workdir: Working directory. If None, a temporary directory is created.
            manager: |TaskManager| object. If None, the object is initialized from the configuration file.

        Returns:
            |PhbstFile| with the phonon band structure.
            |PhdosFile| with the the phonon DOS.
        """
        if ngqpt is None: ngqpt = self.guessed_ngqpt

        if lo_to_splitting == "automatic":
            lo_to_splitting = self.has_lo_to_data() and dipdip != 0

        if lo_to_splitting and not self.has_lo_to_data():
            cprint("lo_to_splitting is True but Eps_inf and Becs are not available in DDB: %s" % self.filepath, "yellow")

        inp = AnaddbInput.phbands_and_dos(
            self.structure, ngqpt=ngqpt, ndivsm=ndivsm, line_density=line_density,
            nqsmall=nqsmall, qppa=qppa, q1shft=(0, 0, 0), qptbounds=qptbounds,
            asr=asr, chneut=chneut, dipdip=dipdip, dos_method=dos_method, lo_to_splitting=lo_to_splitting,
            anaddb_kwargs=anaddb_kwargs, spell_check=spell_check)

        task = self._run_anaddb_task(inp, mpi_procs, workdir, manager, verbose)

        # Open file and add metadata to phbands from DDB
        # TODO: in principle phbands.add_params?
        phbst_file = task.open_phbst()
        self._add_params(phbst_file.phbands)
        if lo_to_splitting:
            phbst_file.phbands.read_non_anal_from_file(os.path.join(task.workdir, "anaddb.nc"))

        phdos_file = None if inp["prtdos"] == 0 else task.open_phdos()
        #if phdos_file is not None: self._add_params(phdos_file.phdos)

        return phbst_file, phdos_file

    def get_coarse(self, ngqpt_coarse, filepath=None):
        """
        Get a version of this file on a coarse mesh

        Args:
            ngqpt_coarse: list of ngqpt indexes that must be a sub-mesh of the original ngqpt
            filepath: Filename for coarse DDB. If None, temporary filename is used.

        Return: DdbFile on coarse mesh.
        """
        # Check if ngqpt is a sub-mesh of ngqpt
        ngqpt_fine = self.guessed_ngqpt
        if any([a % b for a, b in zip(ngqpt_fine, ngqpt_coarse)]):
            raise ValueError('Coarse q-mesh is not a sub-mesh of the current q-mesh')

        # Get the points in the fine mesh
        fine_qpoints = [q.frac_coords for q in self.qpoints]

        # Generate the points of the coarse mesh
        map_fine_to_coarse = []
        nx,ny,nz = ngqpt_coarse
        for i,j,k in itertools.product(range(-int(nx/2), int(nx/2) + 1),
                                       range(-int(ny/2), int(ny/2) + 1),
                                       range(-int(nz/2), int(nz/2) + 1)):
            coarse_qpt = np.array([i, j, k]) / np.array(ngqpt_coarse)
            for n,fine_qpt in enumerate(fine_qpoints):
                if np.allclose(coarse_qpt, fine_qpt):
                    map_fine_to_coarse.append(n)

        # Write the file with a subset of q-points
        if filepath is None:
            import tempfile
            _, filepath = tempfile.mkstemp(suffix="_DDB", text=True)

        self.write(filepath, filter_blocks=map_fine_to_coarse)
        return self.__class__(filepath)

    def anacompare_asr(self, asr_list=(0, 2), chneut_list=(1,), dipdip=1, lo_to_splitting="automatic",
                       nqsmall=10, ndivsm=20, dos_method="tetra", ngqpt=None,
                       verbose=0, mpi_procs=1):
        """
        Invoke anaddb to compute the phonon band structure and the phonon DOS with different
        values of the ``asr`` input variable (acoustic sum rule treatment).
        Build and return |PhononBandsPlotter| object.

        Args:
            asr_list: List of ``asr`` values to test.
            chneut_list: List of ``chneut`` values to test (used by anaddb only if dipdip == 1).
            dipdip: 1 to activate treatment of dipole-dipole interaction (requires BECS and dielectric tensor).
            lo_to_splitting: Allowed values are [True, False, "automatic"]. Defaults to "automatic"
                If True the LO-TO splitting will be calculated if qpoint == Gamma and the non_anal_directions
                non_anal_phfreqs attributes will be addeded to the phonon band structure.
                "automatic" activates LO-TO if the DDB file contains the dielectric tensor and Born effective charges.
            nqsmall: Defines the q-mesh for the phonon DOS in terms of
                the number of divisions to be used to sample the smallest reciprocal lattice vector.
                0 to disable DOS computation.
            ndivsm: Number of division used for the smallest segment of the q-path
            dos_method: Technique for DOS computation in  Possible choices: "tetra", "gaussian" or "gaussian:0.001 eV".
                In the later case, the value 0.001 eV is used as gaussian broadening
            ngqpt: Number of divisions for the ab-initio q-mesh in the DDB file. Auto-detected if None (default)
            verbose: Verbosity level.
            mpi_procs: Number of MPI processes used by anaddb.

        Return:
            |PhononBandsPlotter| object.

            Client code can use ``plotter.combiplot()`` or ``plotter.gridplot()``
            to visualize the results.
        """
        phbands_plotter = PhononBandsPlotter()

        for asr, chneut in itertools.product(asr_list, chneut_list):
            phbst_file, phdos_file = self.anaget_phbst_and_phdos_files(
                nqsmall=nqsmall, ndivsm=ndivsm, asr=asr, chneut=chneut, dipdip=dipdip, dos_method=dos_method,
                lo_to_splitting=lo_to_splitting, ngqpt=ngqpt, qptbounds=None,
                anaddb_kwargs=None, verbose=verbose, mpi_procs=mpi_procs, workdir=None, manager=None)

            label = "asr: %d, dipdip: %d, chneut: %d" % (asr, dipdip, chneut)
            if phdos_file is not None:
                phbands_plotter.add_phbands(label, phbst_file.phbands, phdos=phdos_file.phdos)
                phdos_file.close()
            else:
                phbands_plotter.add_phbands(label, phbst_file.phbands)
            phbst_file.close()

        return phbands_plotter

    def anacompare_dipdip(self, chneut_list=(1,), asr=2, lo_to_splitting="automatic",
                          nqsmall=10, ndivsm=20, dos_method="tetra", ngqpt=None,
                          verbose=0, mpi_procs=1):
        """
        Invoke anaddb to compute the phonon band structure and the phonon DOS with different
        values of the ``asr`` input variable (acoustic sum rule treatment).
        Build and return |PhononDosPlotter| object.

        Args:
            chneut_list: List of ``chneut`` values to test (used for dipdip == 1).
            lo_to_splitting: Allowed values are [True, False, "automatic"]. Defaults to "automatic"
                If True the LO-TO splitting will be calculated if qpoint == Gamma and the non_anal_directions
                non_anal_phfreqs attributes will be addeded to the phonon band structure.
                "automatic" activates LO-TO if the DDB file contains the dielectric tensor and Born effective charges.
            nqsmall: Defines the q-mesh for the phonon DOS in terms of
                the number of divisions to be used to sample the smallest reciprocal lattice vector.
                0 to disable DOS computation.
            ndivsm: Number of division used for the smallest segment of the q-path
            dos_method: Technique for DOS computation in  Possible choices: "tetra", "gaussian" or "gaussian:0.001 eV".
                In the later case, the value 0.001 eV is used as gaussian broadening
            ngqpt: Number of divisions for the ab-initio q-mesh in the DDB file. Auto-detected if None (default)
            verbose: Verbosity level.
            mpi_procs: Number of MPI processes used by anaddb.

        Return:
            |PhononDosPlotter| object.

            Client code can use ``plotter.combiplot()`` or ``plotter.gridplot()``
            to visualize the results.
        """
        phbands_plotter = PhononBandsPlotter()

        for dipdip in (0, 1):
            my_chneut_list = chneut_list if dipdip != 0 else [0]
            for chneut in my_chneut_list:
                phbst_file, phdos_file = self.anaget_phbst_and_phdos_files(
                    nqsmall=nqsmall, ndivsm=ndivsm, asr=asr, chneut=chneut, dipdip=dipdip, dos_method=dos_method,
                    lo_to_splitting=lo_to_splitting, ngqpt=ngqpt, qptbounds=None,
                    anaddb_kwargs=None, verbose=verbose, mpi_procs=mpi_procs, workdir=None, manager=None)

                label = "asr: %d, dipdip: %d, chneut: %d" % (asr, dipdip, chneut)
                if phdos_file is not None:
                    phbands_plotter.add_phbands(label, phbst_file.phbands, phdos=phdos_file.phdos)
                    phdos_file.close()
                else:
                    phbands_plotter.add_phbands(label, phbst_file.phbands)
                phbst_file.close()

        return phbands_plotter

    def anacompare_phdos(self, nqsmalls, asr=2, chneut=1, dipdip=1, dos_method="tetra", ngqpt=None,
                         verbose=0, num_cpus=1, stream=sys.stdout):
        """
        Invoke Anaddb to compute Phonon DOS with different q-meshes. The ab-initio dynamical matrix
        reported in the DDB_ file will be Fourier-interpolated on the list of q-meshes specified
        by ``nqsmalls``. Useful to perform convergence studies.

        Args:
            nqsmalls: List of integers defining the q-mesh for the DOS. Each integer gives
                the number of divisions to be used to sample the smallest reciprocal lattice vector.
            asr, chneut, dipdip: Anaddb input variable. See official documentation.
            dos_method: Technique for DOS computation in  Possible choices: "tetra", "gaussian" or "gaussian:0.001 eV".
                In the later case, the value 0.001 eV is used as gaussian broadening
            ngqpt: Number of divisions for the ab-initio q-mesh in the DDB file. Auto-detected if None (default)
            verbose: Verbosity level.
            num_cpus: Number of CPUs (threads) used to parallellize the calculation of the DOSes. Autodetected if None.
            stream: File-like object used for printing.

        Return:
            ``namedtuple`` with the following attributes::

                    phdoses: List of |PhononDos| objects
                    plotter: |PhononDosPlotter| object.
                        Client code can use ``plotter.gridplot()`` to visualize the results.
        """
        num_cpus = get_ncpus() // 2 if num_cpus is None else num_cpus
        if num_cpus <= 0: num_cpus = 1
        num_cpus = min(num_cpus, len(nqsmalls))

        def do_work(nqsmall):
            phbst_file, phdos_file = self.anaget_phbst_and_phdos_files(
                nqsmall=nqsmall, ndivsm=1, asr=asr, chneut=chneut, dipdip=dipdip, dos_method=dos_method, ngqpt=ngqpt)
            phdos = phdos_file.phdos
            phbst_file.close()
            phdos_file.close()
            return phdos

        if num_cpus == 1:
            # Sequential version
            phdoses = [do_work(nqs) for nqs in nqsmalls]

        else:
            # Threads
            if verbose:
                print("Computing %d phonon DOS with %d threads" % (len(nqsmalls), num_cpus) )
            phdoses = [None] * len(nqsmalls)

            def worker():
                while True:
                    nqsm, phdos_index = q.get()
                    phdos = do_work(nqsm)
                    phdoses[phdos_index] = phdos
                    q.task_done()

            from threading import Thread
            try:
                from Queue import Queue # py2k
            except ImportError:
                from queue import Queue # py3k

            q = Queue()
            for i in range(num_cpus):
                t = Thread(target=worker)
                t.daemon = True
                t.start()

            for i, nqsmall in enumerate(nqsmalls):
                q.put((nqsmall, i))

            # block until all tasks are done
            q.join()

        # Compute relative difference wrt last phonon DOS. Be careful because the DOSes may be defined
        # on different frequency meshes ==> spline on the mesh of the last DOS.
        last_mesh, converged = phdoses[-1].mesh, False
        for i, phdos in enumerate(phdoses[:-1]):
            splined_dos = phdos.spline_on_mesh(last_mesh)
            abs_diff = (splined_dos - phdoses[-1]).abs()
            if verbose:
                print(" Delta(Phdos[%d] - Phdos[%d]) / Phdos[%d]: %f" %
                    (i, len(phdoses)-1, len(phdoses)-1, abs_diff.integral().values[-1]), file=stream)

        # Fill the plotter.
        plotter = PhononDosPlotter()
        for nqsmall, phdos in zip(nqsmalls, phdoses):
            plotter.add_phdos(label="nqsmall %d" % nqsmall, phdos=phdos)

        return dict2namedtuple(phdoses=phdoses, plotter=plotter)

    def anaget_epsinf_and_becs(self, chneut=1, mpi_procs=1, workdir=None, manager=None, verbose=0):
        """
        Call anaddb to compute the macroscopic electronic dielectric tensor (e_inf)
        in Cartesian coordinates and the Born effective charges.

        Args:
            chneut: Anaddb input variable. See official documentation.
            manager: |TaskManager| object. If None, the object is initialized from the configuration file
            mpi_procs: Number of MPI processes to use.
            verbose: verbosity level. Set it to a value > 0 to get more information

        Return: ``namedtuple`` with the following attributes::
            epsinf: |DielectricTensor| object.
            becs: Becs objects.
        """
        if not self.has_lo_to_data():
            cprint("Dielectric tensor and Becs are not available in DDB: %s" % self.filepath, "yellow")

        inp = AnaddbInput(self.structure, anaddb_kwargs={"chneut": chneut})

        task = self._run_anaddb_task(inp, mpi_procs, workdir, manager, verbose)

        # Read data from the netcdf output file produced by anaddb.
        with ETSF_Reader(os.path.join(task.workdir, "anaddb.nc")) as r:
            epsinf = DielectricTensor(r.read_value("emacro_cart").T.copy())
            structure = r.read_structure()
            becs = Becs(r.read_value("becs_cart"), structure, chneut=inp["chneut"], order="f")
            return dict2namedtuple(epsinf=epsinf, becs=becs)

    @deprecated(message="anaget_emacro_and_becs is deprecated and will be removed in abipy 0.8, use anaget_epsinf_and_becs")
    def anaget_emacro_and_becs(self, **kwargs):
        r = self.anaget_epsinf_and_becs(**kwargs)
        return r.epsinf, r.becs

    def anaget_ifc(self, ifcout=None, asr=2, chneut=1, dipdip=1, ngqpt=None,
                   mpi_procs=1, workdir=None, manager=None, verbose=0, anaddb_kwargs=None):
        """
        Execute anaddb to compute the interatomic forces.

        Args:
            ifcout: Number of neighbouring atoms for which the ifc's will be output. If None all the atoms in the big box.
            asr, chneut, dipdip: Anaddb input variable. See official documentation.
            ngqpt: Number of divisions for the q-mesh in the DDB file. Auto-detected if None (default)
            mpi_procs: Number of MPI processes to use.
            workdir: Working directory. If None, a temporary directory is created.
            manager: |TaskManager| object. If None, the object is initialized from the configuration file
            verbose: verbosity level. Set it to a value > 0 to get more information
            anaddb_kwargs: additional kwargs for anaddb

        Returns:
            :class:`InteratomicForceConstants` with the calculated ifc.
        """
        if ngqpt is None: ngqpt = self.guessed_ngqpt

        inp = AnaddbInput.ifc(self.structure, ngqpt=ngqpt, ifcout=ifcout, q1shft=(0, 0, 0), asr=asr, chneut=chneut,
                              dipdip=dipdip, anaddb_kwargs=anaddb_kwargs)

        task = self._run_anaddb_task(inp, mpi_procs, workdir, manager, verbose)

        return InteratomicForceConstants.from_file(os.path.join(task.workdir, 'anaddb.nc'))

    def anaget_dielectric_tensor_generator(self, asr=2, chneut=1, dipdip=1, workdir=None, mpi_procs=1,
                                           manager=None, verbose=0, anaddb_kwargs=None):
        """
        Execute anaddb to extract the quantities necessary to create a |DielectricTensorGenerator|.
        Requires phonon perturbations at Gamma and static electric field perturbations.

        Args:
            asr, chneut, dipdip: Anaddb input variable. See official documentation.
            workdir: Working directory. If None, a temporary directory is created.
            mpi_procs: Number of MPI processes to use.
            manager: |TaskManager| object. If None, the object is initialized from the configuration file
            verbose: verbosity level. Set it to a value > 0 to get more information
            anaddb_kwargs: additional kwargs for anaddb

        Return: |DielectricTensorGenerator| object.
        """
        # Check if gamma is in the DDB.
        try:
            self.qindex((0,0,0))
        except:
            raise ValueError("Gamma point not in %s.\nddb.qpoints:\n%s" % (self.filepath, self.qpoints))

        inp = AnaddbInput.modes_at_qpoint(self.structure, (0, 0, 0), asr=asr, chneut=chneut, dipdip=dipdip,
                                          lo_to_splitting=False, anaddb_kwargs=anaddb_kwargs)

        if anaddb_kwargs is None or 'dieflag' not in anaddb_kwargs:
            inp['dieflag'] = 1

        task = self._run_anaddb_task(inp, mpi_procs, workdir, manager, verbose)

        return DielectricTensorGenerator.from_files(os.path.join(task.workdir, "run.abo_PHBST.nc"),
                                                    os.path.join(task.workdir, "anaddb.nc"))

    def anaget_elastic(self, relaxed_ion="automatic", piezo="automatic",
                        dde=False, stress_correction=False, asr=2, chneut=1,
                        mpi_procs=1, workdir=None, manager=None, verbose=0, retpath=False):
        """
        Call anaddb to compute elastic and piezoelectric tensors. Require DDB with strain terms.

	By default, this method sets the anaddb input variables automatically
	by looking at the 2nd-order derivatives available in the DDB file.
	This behaviour can be changed by setting explicitly the value of:
        `relaxed_ion` and `piezo`.

        Args:
            relaxed_ion: Activate computation of relaxed-ion tensors.
		Allowed values are [True, False, "automatic"]. Defaults to "automatic".
                In "automatic" mode, relaxed-ion tensors are automatically computed if
                internal strain terms and phonons at Gamma are present in the DDB.
            piezo: Activate computation of piezoelectric tensors.
		Allowed values are [True, False, "automatic"]. Defaults to "automatic".
                In "automatic" mode, piezoelectric tensors are automatically computed if
                piezoelectric terms are present in the DDB.
                NB: relaxed-ion piezoelectric requires the activation of `relaxed_ion`.
            dde: if True, dielectric tensors will be calculated.
	    stress_correction: Calculate the relaxed ion elastic tensors, considering
                the stress left inside cell. The DDB must contain the stress tensor.
            asr: Anaddb input variable. See official documentation.
            chneut: Anaddb input variable. See official documentation.
            mpi_procs: Number of MPI processes to use.
            workdir: Working directory. If None, a temporary directory is created.
            manager: |TaskManager| object. If None, the object is initialized from the configuration file
            verbose: verbosity level. Set it to a value > 0 to get more information
            retpath: True to return path to anaddb.nc file.

        Return:
            |ElasticData| object if `retpath` is None else absolute path to anaddb.nc file.
        """
        if not self.has_strain_terms(): # DOH!
            cprint("Strain perturbations are not available in DDB: %s" % self.filepath, "yellow")

        if relaxed_ion == "automatic":
            relaxed_ion = self.has_internalstrain_terms() and self.has_at_least_one_atomic_perturbation(qpt=(0, 0, 0))

        if relaxed_ion:
            if not self.has_at_least_one_atomic_perturbation(qpt=(0, 0, 0)):
                cprint("Requiring `relaxed_ion` but no atomic term available in DDB: %s" % self.filepath, "yellow")
            if not self.has_internalstrain_terms():
                cprint("Requiring `internal_strain` but no internal strain term in DDB: %s" % self.filepath, "yellow")

        if piezo == "automatic":
            piezo = self.has_piezoelectric_terms()

        if piezo and not self.has_piezoelectric_terms():
            cprint("Requiring `piezo` but no piezoelectric term available in DDB: %s" % self.filepath, "yellow")

	# FIXME This is problematic so don't use automatic as default
        #select = "all"
        select = "at_least_one_diagoterm"
        if dde == "automatic":
            dde = self.has_epsinf_terms(select=select)

        if dde and not self.has_epsinf_terms(select=select):
            cprint("Requiring `dde` but dielectric tensor not available in DDB: %s" % self.filepath, "yellow")

        if stress_correction == "automatic":
            stress_correction = self.cart_stress_tensor is not None

        if stress_correction and self.cart_stress_tensor is None:
            cprint("Requiring `stress_correction` but stress not available in DDB: %s" % self.filepath, "yellow")

        inp = AnaddbInput.dfpt(self.structure, strain=True, relaxed_ion=relaxed_ion,
		               dde=dde, piezo=piezo, stress_correction=stress_correction, dte=False,
                               asr=asr, chneut=chneut)

        task = self._run_anaddb_task(inp, mpi_procs, workdir, manager, verbose)

        # Read data from the netcdf output file produced by anaddb.
        path = os.path.join(task.workdir, "anaddb.nc")
        return ElasticData.from_file(path) if not retpath else path

    def _run_anaddb_task(self, anaddb_input, mpi_procs, workdir, manager, verbose):
        """
        Execute an |AnaddbInput| via the shell. Return AnaddbTask.
        """
        task = AnaddbTask.temp_shell_task(anaddb_input, ddb_node=self.filepath,
                mpi_procs=mpi_procs, workdir=workdir, manager=manager)

        if verbose:
            print("ANADDB INPUT:\n", anaddb_input)
            print("workdir:", task.workdir)

        # Run the task here.
        task.start_and_wait(autoparal=False)

        report = task.get_event_report()
        if not report.run_completed:
            raise self.AnaddbError(task=task, report=report)

        return task

    def write(self, filepath, filter_blocks=None):
        """
        Writes the DDB file in filepath. Requires the blocks data.
        Only the information stored in self.header.lines and in self.blocks will be used to produce the file
        """
        lines = list(self.header.lines)

        if filter_blocks is None:
            blocks = self.blocks
        else:
            blocks = [self.blocks[i] for i in filter_blocks]

        lines.append(" **** Database of total energy derivatives ****")
        lines.append(" Number of data blocks={0:5}".format(len(blocks)))
        lines.append(" ")

        for b in blocks:
            lines.extend(b["data"])
            lines.append(" ")

        lines.append(" List of bloks and their characteristics")
        lines.append(" ")

        for b in blocks:
            lines.extend(b["data"][:2])
            lines.append(" ")

        with open(filepath, "wt") as f:
            f.write("\n".join(lines))

    def get_block_for_qpoint(self, qpt):
        """
        Extracts the block data for the selected qpoint.
        Returns a list of lines containing the block information
        """
        if hasattr(qpt, "frac_coords"): qpt = qpt.frac_coords

        for b in self.blocks:
            if b['qpt'] is not None and np.allclose(b['qpt'], qpt):
                return b["data"]

    def replace_block_for_qpoint(self, qpt, data):
        """
        Change the block data for the selected qpoint. Object is modified in-place.
        Data should be a list of strings representing the whole data block.
        Note that the DDB file should be written and reopened in order to call methods that run anaddb.

        Return:
            True if qpt has been found and data has been replaced.
        """
        if hasattr(qpt, "frac_coords"): qpt = qpt.frac_coords

        for b in self.blocks:
            if b['qpt'] is not None and np.allclose(b['qpt'], qpt):
                b["data"] = data
                return True

        return False

    def write_notebook(self, nbpath=None):
        """
        Write an jupyter_ notebook to nbpath. If ``nbpath`` is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ddb = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("units = 'eV'\nprint(ddb)"),
            nbv.new_code_cell("# display(ddb.header)"),
            nbv.new_markdown_cell("## Invoke `anaddb` to compute bands and dos"),
            nbv.new_code_cell("""\
bstfile, phdosfile =  ddb.anaget_phbst_and_phdos_files(nqsmall=10, ndivsm=20,
    asr=2, chneut=1, dipdip=0, lo_to_splitting="automatic",
    dos_method="tetra", ngqpt=None, qptbounds=None, verbose=0, anaddb_kwargs=None)

phbands, phdos = bstfile.phbands, phdosfile.phdos"""),
            nbv.new_markdown_cell("## q-point path"),
            nbv.new_code_cell("phbands.qpoints.plot();"),
            nbv.new_markdown_cell("## Phonon bands with DOS"),
            nbv.new_code_cell("phbands.plot_with_phdos(phdos, units=units);"),
            nbv.new_markdown_cell("## Phonon fatbands with DOS"),
            nbv.new_code_cell("phbands.plot_fatbands(phdos_file=phdosfile, units=units);"),
            nbv.new_markdown_cell("## Distribution of phonon frequencies wrt mode index"),
            nbv.new_code_cell("phbands.boxplot(units=units);"),
            nbv.new_markdown_cell("## Phonon band structure with different color for each line"),
            nbv.new_code_cell("phbands.plot_colored_matched(units=units);"),
            nbv.new_markdown_cell("## Type-projected phonon DOS."),
            nbv.new_code_cell("phdosfile.plot_pjdos_type(units=units);"),
            nbv.new_markdown_cell("## Type-projected phonon DOS decomposed along the three reduced directions"),
            nbv.new_code_cell("phdosfile.plot_pjdos_redirs_type(units=units);"),
            nbv.new_code_cell("#phdosfile.plot_pjdos_redirs_site(units=units);"),
            nbv.new_markdown_cell("## Thermodinamic properties within the harmonic approximation"),
            nbv.new_code_cell("phdosfile.phdos.plot_harmonic_thermo(tstart=5, tstop=300);"),

            nbv.new_markdown_cell("## Macroscopic dielectric tensor and Born effective charges"),
            nbv.new_code_cell("""\
if False:
    einf, becs = ddb.anaget_epsinf_and_becs()
    print(einf)
    print(becs)"""),

            nbv.new_markdown_cell("## Call `anaddb` to compute phonons and DOS with/without ASR"),
            nbv.new_code_cell("""\
#asr_plotter = ddb.anacompare_asr(asr_list=(0, 2), nqsmall=0, ndivsm=10)
#asr_plotter.gridplot();
"""),

            nbv.new_markdown_cell("## Call `anaddb` to compute phonon DOS with different BZ samplings"),
            nbv.new_code_cell("""\
c = None
if False:
    c = ddb.anacompare_phdos(nqsmalls=[20, 30], asr=2, chneut=1, dipdip=1,
            dos_method="tetra", ngqpt=None, num_cpus=1)"""),

            nbv.new_code_cell("""\
if c is not None:
    c.plotter.ipw_select_plot()
    #for phdos in c.phdoses"""),

            nbv.new_markdown_cell("## Analysis of the IFCs in real-space"),
            nbv.new_code_cell("""\
ifc = None
if False:
    ifc = ddb.anaget_ifc(ifcout=None, asr=2, chneut=1, dipdip=1, ngqpt=None, verbose=0, anaddb_kwargs=None)"""),

            nbv.new_code_cell("""\
if ifc is not None:
    ifc.plot_longitudinal_ifc(atom_indices=None, atom_element=None, neighbour_element=None, min_dist=None,
                                    max_dist=None, ax=None);"""),
            nbv.new_code_cell("""\
if ifc is not None:
    ifc.plot_longitudinal_ifc_short_range(atom_indices=None, atom_element=None, neighbour_element=None);"""),
            nbv.new_code_cell("""\
if ifc is not None:
    ifc.plot_longitudinal_ifc_ewald(atom_indices=None, atom_element=None, neighbour_element=None,
                                    min_dist=None, max_dist=None, ax=None);"""),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class Becs(Has_Structure):
    """
    This object stores the Born effective charges and provides simple tools for data analysis.
    """

    def __init__(self, becs_arr, structure, chneut, order="c"):
        """
        Args:
            becs_arr: [3, 3, natom] array with the Born effective charges in Cartesian coordinates.
            structure: |Structure| object.
            chneut: Option used for the treatment of the Charge Neutrality requirement
                for the effective charges (anaddb input variable)
            order: "f" if becs_arr is in Fortran order.
        """
        assert len(becs_arr) == len(structure)
        self._structure = structure
        self.chneut = chneut

        # Values is a numpy array while zstars is a list of Tensor objects.
        self.values = np.empty((len(structure), 3, 3))
        for i, bec in enumerate(becs_arr):
            mat = becs_arr[i]
            if order.lower() == "f": mat = mat.T.copy()
            self.values[i] = mat

        self.zstars = [ZstarTensor(mat) for mat in self.values]

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def __repr__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app("Born effective charges in Cartesian coordinates (Voigt notation)")
        app(self.get_voigt_dataframe().to_string())
        app("")

        if verbose:
            app("Born effective charges (full tensor)")
            for site, bec in zip(self.structure, self.values):
                app("Z* at site: %s" % repr(site))
                app(str(bec))
                app("")

        # Add info on the bec sum rule.
        app("Born effective charge neutrality sum-rule with chneut: %d\n" % self.chneut)
        app(str(self.sumrule))

        return "\n".join(lines)

    @property
    def sumrule(self):
        """[3, 3] matrix with Born effective charge neutrality sum-rule."""
        return self.values.sum(axis=0)

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        return self.get_voigt_dataframe()._repr_html_()

    def get_voigt_dataframe(self, view="inequivalent", tol=1e-3, select_symbols=None, decimals=5, verbose=0):
        """
        Return |pandas-DataFrame| with Voigt indices as columns and natom rows.

        Args:
            view: "inequivalent" to show only inequivalent atoms. "all" for all sites.
            tol: Entries are set to zero below this value
            select_symbols: String or list of strings with chemical symbols.
                Used to select only atoms of this type.
            decimals: Number of decimal places to round to.
                If decimals is negative, it specifies the number of positions to the left of the decimal point.
            verbose: Verbosity level.
        """
        aview = self._get_atomview(view, select_symbols=select_symbols, verbose=verbose)

        columns = ["xx", "yy", "zz", "yz", "xz", "xy"]
        rows = []
        for (iatom, wlabel) in zip(aview.iatom_list, aview.wyck_labels):
            site = self.structure[iatom]
            zstar = self.zstars[iatom]
            d = OrderedDict()
            d["element"] = site.specie.symbol
            d["site_index"] = iatom
            d["frac_coords"] = np.round(site.frac_coords, decimals=decimals)
            d["cart_coords"] = np.round(site.coords, decimals=decimals)
            d["wyckoff"] = wlabel
            zstar = zstar.zeroed(tol=tol)
            for k, v in zip(columns, zstar.voigt):
                d[k] = v
            if verbose:
                d["determinant"] = np.linalg.det(zstar)
                d["iso"] = zstar.trace() / 3
            rows.append(d)

        return pd.DataFrame(rows, columns=list(rows[0].keys()) if rows else None)

    def check_site_symmetries(self, verbose=0):
        from abipy.core.wyckoff import SiteSymmetries
        ss = SiteSymmetries(self.structure)
        return ss.check_site_symmetries(self.values, verbose=verbose)


class DielectricTensorGenerator(Has_Structure):
    """
    Object used to generate frequency dependent dielectric tensors as obtained
    from DFPT calculations. The values are calculated on the fly
    based on the phonon frequencies at gamma and oscillator strengths.
    The first three frequencies would be considered as acoustic modes and
    ignored in the calculation. No checks would be performed.

    See the definitions Eq.(53-54) in :cite:`Gonze1997` PRB55, 10355 (1997).
    """

    @classmethod
    def from_files(cls, phbst_filepath, anaddbnc_filepath):
        """
        Generates the object from the files that contain the phonon frequencies, oscillator strength and
        static dielectric tensor, i.e. the PHBST.nc and anaddb.nc netcdf files, respectively.
        """
        with ETSF_Reader(phbst_filepath) as reader:
            qpts = reader.read_value("qpoints")
            full_phfreqs = reader.read_value("phfreqs")

        for i, q in enumerate(qpts):
            if np.array_equal(q, [0, 0, 0]):
                phfreqs = full_phfreqs[i].copy()
                break
        else:
            raise ValueError('The PHBST does not contain frequencies at gamma')

        with ETSF_Reader(anaddbnc_filepath) as reader:
            epsinf = DielectricTensor(reader.read_value("emacro_cart").T.copy())
            eps0 = DielectricTensor(reader.read_value("emacro_cart_rlx").T.copy())
            try:
                oscillator_strength = reader.read_value("oscillator_strength", cmode="c")
                oscillator_strength = oscillator_strength.transpose((0, 2, 1)).copy()
            except Exception as exc:
                import traceback
                msg = traceback.format_exc()
                msg += ("Error while trying to read from file.\n"
                        "Verify that dieflag == 1, 3 or 4 in anaddb\n")
                raise ValueError(msg)

            structure = reader.read_structure()

        return cls(phfreqs, oscillator_strength, eps0, epsinf, structure)

    @classmethod
    def from_objects(cls, phbands, anaddbnc):
        """
        Generates the object from the objects |PhononBands| object and AnaddbNcFile
        """
        gamma_index = phbands.qindex([0, 0, 0])

        phfreqs = phbands.phfreqs[gamma_index]

        epsinf = anaddbnc.epsinf
        eps0 = anaddbnc.eps0
        oscillator_strength = anaddbnc.oscillator_strength

        return cls(phfreqs, oscillator_strength, eps0, epsinf, anaddbnc.structure)

    def __init__(self, phfreqs, oscillator_strength, eps0, epsinf, structure):
        """
        Args:
             phfreqs: numpy array containing the 3 * num_atoms phonon frequencies at gamma
             oscillator_strength: complex numpy array with shape [number of phonon modes, 3, 3] in atomic units
             eps0: numpy array containing the e0 dielectric tensor without frequency dependence
             epsinf: numpy array with the electronic dielectric tensor (einf) without frequency dependence
             structure: |Structure| object.
        """
        self.phfreqs = phfreqs
        self.oscillator_strength = oscillator_strength
        self.eps0 = eps0
        self.epsinf = epsinf
        self._structure = structure

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        lines = []
        app = lines.append
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        app(marquee("Oscillator strength", mark="="))
        tol = 1e-8
        app("Real part in Cartesian coordinates. a.u. units; 1 a.u. = 253.2638413 m3/s2. Set to zero below %.2e." % tol)
        app(self.get_oscillator_dataframe(reim="re", tol=tol).to_string())
        if verbose:
            app("")
            app("Imaginary part in a.u.; 1 a.u. = 253.2638413 m3/s2. Set to zero below %.2e." % tol)
            app(self.get_oscillator_dataframe(reim="im", tol=tol).to_string())
            app("")
            app("Trace of oscillator strength, for each phonon mode:")
            traces = [o.trace() for o in self.oscillator_strength]
            app(str(traces))
        app("")

        tol = 1e-3
        app(marquee("Dielectric Tensors", mark="="))
        app("Electronic dielectric tensor (eps_inf) in Cartesian coordinates. Set to zero below %.2e." % tol)
        app(self.epsinf.get_dataframe(tol=tol).to_string())
        app("")
        app("Zero-frequency dielectric tensor (eps_zero) in Cartesian coordinates. Set to zero below %.2e." % tol)
        app(self.eps0.get_dataframe(tol=tol).to_string())

        return "\n".join(lines)

    def get_oscillator_dataframe(self, reim="all", tol=1e-8):
        """
        Return |pandas-Dataframe| with oscillator matrix elements.

        Args:
            reim: "re" for real part, "im" for imaginary part, "all" for both.
            tol: Entries are set to zero below this value
        """
        dmap = dict(xx=(0, 0), yy=(1, 1), zz=(2, 2), yz=(1, 2), xz=(0, 2), xy=(0, 1))
        #decimals = int(abs(np.rint(np.log10(tol))))
        # 1 a.u. = 253.2638413 m3/s2.
        # TODO: Use SI?
        #fact = 253.2638413

        rows, index = [], []
        for nu in range(3 * len(self.structure)):
            d = {k: data_from_cplx_mode(reim, self.oscillator_strength[nu][t], tol=tol) for k, t in dmap.items()}
            #d = {k: np.around(v * fact, decimals=decimals) for k, v in d.items()}
            rows.append(d)
            index.append(nu)

        df = pd.DataFrame(rows, index=index, columns=list(rows[0].keys()))
        df.index.name = "mode"
        return df

    def tensor_at_frequency(self, w, gamma_ev=1e-4, units='eV'):
        """
        Returns a |DielectricTensor| object representing the dielectric tensor
        in atomic units at the specified frequency w. Eq.(53-54) in PRB55, 10355 (1997).

        Args:
            w: Frequency in eV
            gamma_ev: Phonon damping factor in eV (full width). Poles are shifted by phfreq * gamma_ev.
                Accept scalar or [nfreq] array.
            units: string specifying the units used for ph frequencies.  Possible values in
            ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
        """
        w =  w / phfactor_ev2units(units)

        # Note that the acoustic modes are not included: their oscillator strength should be exactly zero
        # Also, only the real part of the oscillators is taken into account:
        # the possible imaginary parts of degenerate modes will cancel.
        if duck.is_listlike(gamma_ev):
            gammas = np.asarray(gamma_ev)
            assert len(gammas) == len(phfreqs)
        else:
            gammas = np.ones(len(self.phfreqs)) * float(gamma_ev)

        t = np.zeros((3, 3),dtype=complex)
        for i in range(3, len(self.phfreqs)):
            g =  gammas[i] * self.phfreqs[i]
            t += self.oscillator_strength[i].real / (self.phfreqs[i]**2 - w**2 - 1j*g)

        vol = self.structure.volume / bohr_to_angstrom ** 3
        t = 4 * np.pi * t / vol / eV_to_Ha**2
        t += self.epsinf

        #t = np.zeros((3, 3))
        #w = w * eV_to_Ha
        #for i in range(3, len(self.phfreqs)):
        #    phw = self.phfreqs[i] * eV_to_Ha
        #    t += self.oscillator_strength[i].real / (phw**2 - w**2)
        #t *= 4 * np.pi / (self.structure.volume * ang_to_bohr ** 3)
        #t += self.epsinf

        return DielectricTensor(t)

    @add_fig_kwargs
    def plot(self, w_min=0, w_max=None, gamma_ev=1e-4, num=500, component='diag', reim="reim", units='eV',
             with_phfreqs=True, ax=None, fontsize=12, **kwargs):
        """
        Plots the selected components of the dielectric tensor as a function of frequency.

        Args:
            w_min: minimum frequency in units `units`.
            w_max: maximum frequency. If None it will be set to the value of the maximum frequency + 5*gamma_ev.
            gamma_ev: Phonon damping factor in eV (full width). Poles are shifted by phfreq * gamma_ev.
                Accept scalar or [nfreq] array.
            num: number of values of the frequencies between w_min and w_max.
            component: determine which components of the tensor will be displayed. Can be a list/tuple of two
                elements, indicating the indices [i, j] of the desired component or a string among:

                * 'diag_av' to plot the average of the components on the diagonal
                * 'diag' to plot the elements on diagonal
                * 'all' to plot all the components in the upper triangle.
                * 'offdiag' to plot the off-diagonal components in the upper triangle.

            reim: a string with "re" will plot the real part, with "im" selects the imaginary part.
            units: string specifying the units used for phonon frequencies. Possible values in
                ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            with_phfreqs: True to show phonon frequencies with dots.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: Legend and label fontsize.

        Return: |matplotlib-Figure|
        """
        if w_max is None:
            w_max = np.max(self.phfreqs) * phfactor_ev2units(units) + gamma_ev * 10

        wmesh = np.linspace(w_min, w_max, num, endpoint=True)
        t = np.zeros((num, 3, 3), dtype=complex)

        for i, w in enumerate(wmesh):
            t[i] = self.tensor_at_frequency(w, units=units, gamma_ev=gamma_ev)

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = 2

        ax.set_xlabel('Frequency {}'.format(phunit_tag(units)))
        ax.set_ylabel(r'$\epsilon(\omega)$')
        ax.grid(True)

        reimfs = []
        if 're' in reim: reimfs.append((np.real, "Re{%s}"))
        if 'im' in reim: reimfs.append((np.imag, "Im{%s}"))

        for reimf, reims in reimfs:
            if isinstance(component, (list, tuple)):
                label = reims % r'$\epsilon_{%d%d}$' % tuple(component)
                ax.plot(wmesh, reimf(t[:,component[0], component[1]]), label=label, **kwargs)
            elif component == 'diag':
                for i in range(3):
                    label = reims % r'$\epsilon_{%d%d}$' % (i, i)
                    ax.plot(wmesh, reimf(t[:, i, i]), label=label, **kwargs)
            elif component in ('all', "offdiag"):
                for i in range(3):
                    for j in range(3):
                        if component == "all" and i > j: continue
                        if component == "offdiag" and i >= j: continue
                        label = reims % r'$\epsilon_{%d%d}$' % (i, j)
                        ax.plot(wmesh, reimf(t[:, i, j]), label=label, **kwargs)
            elif component == 'diag_av':
                label = r'$Average\, %s\epsilon_{ii}$' % reims
                ax.plot(wmesh, np.trace(reimf(t), axis1=1, axis2=2)/3, label=label, **kwargs)
            else:
                raise ValueError('Unkwnown component {}'.format(component))

        # Add points showing phonon energies.
        if with_phfreqs:
            wvals = self.phfreqs[3:] * phfactor_ev2units(units)
            ax.scatter(wvals, np.zeros_like(wvals), s=30, marker="o", c="blue")

        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        return fig

    # To maintain backward compatibility.
    plot_vs_w = plot

    @add_fig_kwargs
    def plot_all(self, **kwargs):
        """
        Plot diagonal and off-diagonal elements of the dielectric tensor as a function of frequency.
        Both real and imag part are show. Accepts all arguments of `plot` method with the exception of:
            `component` and `reim`.

        Returns: |matplotlib-Figure|
        """
        axmat, fig, plt = get_axarray_fig_plt(None, nrows=2, ncols=2,
                                              sharex=True, sharey=False, squeeze=False)
        fontsize = kwargs.pop("fontsize", 8)
        for irow in range(2):
            component = {0: "diag", 1: "offdiag"}[irow]
            for icol in range(2):
                reim = {0: "re", 1: "im"}[icol]
                self.plot(component=component, reim=reim, ax=axmat[irow, icol], fontsize=fontsize, show=False, **kwargs)

        return fig


class DdbRobot(Robot):
    """
    This robot analyzes the results contained in multiple DDB_ files.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: DdbRobot
    """
    EXT = "DDB"

    @classmethod
    def class_handles_filename(cls, filename):
        """Exclude DDB.nc files. Override base class."""
        return filename.endswith("_" + cls.EXT)

    @classmethod
    def from_mpid_list(cls, mpid_list, api_key=None, endpoint=None):
        """
        Build a DdbRobot from list of materials-project ids.

        Args:
            mpid_list: List of Materials Project material_ids (e.g., ["mp-1234", "mp-1245"]).
            api_key (str): A String API key for accessing the MaterialsProject REST interface.
                If None, the code will check if there is a `PMG_MAPI_KEY` in your .pmgrc.yaml.
            endpoint (str): Url of endpoint to access the MaterialsProject REST interface.
                Defaults to the standard Materials Project REST address
        """
        from abipy.core import restapi
        ddb_files = []
        with restapi.get_mprester(api_key=api_key, endpoint=endpoint) as rest:
            for mpid in list_strings(mpid_list):
                try:
                    ddb_string = rest._make_request("/materials/%s/abinit_ddb" % mpid)
                except rest.Error:
                    cprint("Cannot get DDB for mp-id: %s, ignoring error" % mpid, "yellow")
                    continue

                _, tmpfile = tempfile.mkstemp(prefix=mpid, suffix='_DDB')
                ddb_files.append(tmpfile)
                with open(tmpfile, "wt") as fh:
                    fh.write(ddb_string)

        return cls.from_files(ddb_files)

    #def get_qpoints_union(self):
    #    """
    #    Return numpy array with the q-points in reduced coordinates found in the DDB files.
    #    """
    #    qpoints = []
    #    for label, ddb in self.items():
    #        qpoints.extend(q.frac_coords for q in ddb.qpoints if q not in qpoints)

    #    return np.array(qpoints)

    #def get_qpoints_intersection(self):
    #    """Return numpy array with the q-points in reduced coordinates found in the DDB files."""
    #    qpoints = []
    #    for label, ddb in self.items():
    #        qpoints.extend(q.frac_coords for q in ddb.qpoints if q not in qpoints)
    #
    #    return np.array(qpoints)

    # DEBUGGING CODE (do not remove)
    #def find_duplicated_entries(self, std_tol=1e-5, verbose=1):
    #    """
    #    Check for duplicated entries in the list of ddb files

    #    Args:
    #        std_tol: Tolerance on standard deviation
    #        verbose: Verbosity level.

    #    Return: (retcode, results) where results maps qpt --> DataFrame with perts as index.
    #    """
    #    from pprint import pprint

    #    # Build q --> group of dataframes.
    #    from collections import defaultdict
    #    q2dfgroup = defaultdict(list)
    #    for ddb in self.abifiles:
    #        for qpt, df in ddb.computed_dynmat.items():
    #            q2dfgroup[qpt].append(df)

    #    retcode, results = 0, {}
    #    for qpt, dfgroup in q2dfgroup.items():
    #        all_indset = [set(df.index) for df in dfgroup]
    #        # Build union of all dynmat indices with this q
    #        allps = set(all_indset[0]).union(*all_indset)
    #        #allps = set(all_indset[0]).intersection(*all_indset)

    #        index, d_list = [], []
    #        for p in allps:
    #            # Find dataframes with this p
    #            found = [p in index for index in all_indset]
    #            count = found.count(True)
    #            if count == 1: continue
    #            if verbose:
    #                print("Found %s duplicated entries for p: %s" % (count, str(p)))

    #            # Compute stats for this p (complex numbers)
    #            cvalues = []
    #            for f, df in zip(found, dfgroup):
    #                if not f: continue
    #                c = df["cvalue"].loc[[p]]
    #                cvalues.append(c)

    #            cvalues = np.array(cvalues)
    #            norms = np.abs(cvalues)
    #            d = dict(mean=cvalues.mean(), std=cvalues.std(),
    #                     min_norm=norms.min(), max_norm=norms.max(), count=count)

    #            # Print warning if large deviation
    #            #if d["max_norm"]  - d["min_norm"] > 1e-5:
    #            if d["std"] > std_tol:
    #                retcode += 1
    #                cprint("Found std > %s" % std_tol, "red")
    #                pprint(cvalues)
    #            if verbose:
    #                pprint(d)
    #                print(2 * "")

    #            d_list.append(d)
    #            index.append(p)

    #        results[qpt] = pd.DataFrame(d_list, index=index)

    #    return retcode, results

    def get_dataframe_at_qpoint(self, qpoint=None, units="eV", asr=2, chneut=1, dipdip=1,
	    with_geo=True, with_spglib=True, abspath=False, funcs=None):
        """
	Call anaddb to compute the phonon frequencies at a single q-point using the DDB files treated
	by the robot and the given anaddb input arguments. LO-TO splitting is not included.
        Build and return a |pandas-Dataframe| with results

        Args:
            qpoint: Reduced coordinates of the qpoint where phonon modes are computed
            units: string specifying the units used for ph frequencies.  Possible values in
                ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.
            asr, chneut, dipdip: Anaddb input variable. See official documentation.
            with_geo: True if structure info should be added to the dataframe
	    with_spglib: True to compute spglib space group and add it to the DataFrame.
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().
            funcs: Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a |DdbFile| object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.

        Return:
            |pandas-DataFrame|
        """
        # If qpoint is None, all the DDB must contain have the same q-point .
        if qpoint is None:
            if not all(len(ddb.qpoints) == 1 for ddb in self.abifiles):
                raise ValueError("Found more than one q-point in the DDB file. qpoint must be specified")

            qpoint = self[0].qpoints[0]
            if any(np.any(ddb.qpoints[0] != qpoint) for ddb in self.abifiles):
                raise ValueError("All the q-points in the DDB files must be equal")

        rows, row_names = [], []
        for i, (label, ddb) in enumerate(self.items()):
            row_names.append(label)
            d = OrderedDict()
            #d = {aname: getattr(ddb, aname) for aname in attrs}
            #d.update({"qpgap": mdf.get_qpgap(spin, kpoint)})

            # Call anaddb to get the phonon frequencies. Note lo_to_splitting set to False.
            phbands = ddb.anaget_phmodes_at_qpoint(qpoint=qpoint, asr=asr, chneut=chneut,
               dipdip=dipdip, lo_to_splitting=False)
            # [nq, nmodes] array
            freqs = phbands.phfreqs[0, :] * phfactor_ev2units(units)

            d.update({"mode" + str(i): freqs[i] for i in range(len(freqs))})

            # Add convergence parameters
            d.update(ddb.params)

            # Add info on structure.
            if with_geo:
                d.update(phbands.structure.get_dict4pandas(with_spglib=with_spglib))

            # Execute functions.
            if funcs is not None: d.update(self._exec_funcs(funcs, ddb))
            rows.append(d)

        row_names = row_names if not abspath else self._to_relpaths(row_names)
        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    def anaget_phonon_plotters(self, **kwargs):
        r"""
        Invoke anaddb to compute phonon bands and DOS using the arguments passed via \*\*kwargs.
        Collect results and return `namedtuple` with the following attributes:

            phbands_plotter: |PhononBandsPlotter| object.
            phdos_plotter: |PhononDosPlotter| object.
        """
	# TODO: Multiprocessing?
        if "workdir" in kwargs:
            raise ValueError("Cannot specify `workdir` when multiple DDB file are executed.")

        phbands_plotter, phdos_plotter = PhononBandsPlotter(), PhononDosPlotter()

        for label, ddb in self.items():
            # Invoke anaddb to get phonon bands and DOS.
            phbst_file, phdos_file = ddb.anaget_phbst_and_phdos_files(**kwargs)

            # Phonon frequencies with non analytical contributions, if calculated, are saved in anaddb.nc
            # Those results should be fetched from there and added to the phonon bands.
            # lo_to_splitting in ["automatic", True, False] and defaults to automatic.
            if kwargs.get("lo_to_splitting", False):
                anaddb_path = os.path.join(os.path.dirname(phbst_file.filepath), "anaddb.nc")
                phbst_file.phbands.read_non_anal_from_file(anaddb_path)

            phbands_plotter.add_phbands(label, phbst_file, phdos=phdos_file)
            phbst_file.close()
            if phdos_file is not None:
                phdos_plotter.add_phdos(label, phdos=phdos_file.phdos)
                phdos_file.close()

        return dict2namedtuple(phbands_plotter=phbands_plotter, phdos_plotter=phdos_plotter)

    def anacompare_elastic(self, ddb_header_keys=None, with_structure=True, with_spglib=True,
                           with_path=False, manager=None, verbose=0, **kwargs):
        """
        Compute elastic and piezoelectric properties for all DDBs in the robot and build DataFrame.

        Args:
            ddb_header_keys: List of keywords in the header of the DDB file
                whose value will be added to the Dataframe.
            with_structure: True to add structure parameters to the DataFrame.
	    with_spglib: True to compute spglib space group and add it to the DataFrame.
            with_path: True to add DDB path to dataframe
            manager: |TaskManager| object. If None, the object is initialized from the configuration file
            verbose: verbosity level. Set it to a value > 0 to get more information
            kwargs: Keyword arguments passed to `ddb.anaget_elastic`.

        Return: DataFrame and list of ElastData objects.
        """
        ddb_header_keys = [] if ddb_header_keys is None else list_strings(ddb_header_keys)
        df_list, elastdata_list = [], []
        for label, ddb in self.items():
            # Invoke anaddb to compute elastic data.
            edata = ddb.anaget_elastic(verbose=verbose, **kwargs)
            elastdata_list.append(edata)

	    # Build daframe with properties derived from the elastic tensor.
            df = edata.get_elastic_properties_dataframe()

            # Add metadata to the dataframe.
            df["formula"] = ddb.structure.formula
            for k in ddb_header_keys:
                df[k] = ddb.header[k]

            # Add structural parameters to the dataframe.
            if with_structure:
                for skey, svalue in ddb.structure.get_dict4pandas(with_spglib=with_spglib).items():
                    df[skey] = svalue

            # Add path to the DDB file.
            if with_path: df["ddb_path"] = ddb.filepath

            df_list.append(df)

        # Concatenate dataframes.
        return dict2namedtuple(df=pd.concat(df_list, ignore_index=True),
                               elastdata_list=elastdata_list)

    def anacompare_becs(self, ddb_header_keys=None, chneut=1, tol=1e-3, with_path=False, verbose=0):
        """
        Compute Born effective charges for all DDBs in the robot and build DataFrame.
        with Voigt indices as columns + metadata. Useful for convergence studies.

        Args:
            ddb_header_keys: List of keywords in the header of the DDB file
                whose value will be added to the Dataframe.
            chneut: Anaddb input variable. See official documentation.
            tol: Elements below this value are set to zero.
            with_path: True to add DDB path to dataframe
            verbose: verbosity level. Set it to a value > 0 to get more information

        Return: ``namedtuple`` with the following attributes::

            df: DataFrame with Voigt as columns.
            becs_list: list of Becs objects.
        """
        ddb_header_keys = [] if ddb_header_keys is None else list_strings(ddb_header_keys)
        df_list, becs_list = [], []
        for label, ddb in self.items():
            # Invoke anaddb to compute Becs
            _, becs = ddb.anaget_epsinf_and_becs(chneut=chneut, verbose=verbose)
            becs_list.append(becs)
            df = becs.get_voigt_dataframe(tol=tol)

            # Add metadata to the dataframe.
            df["formula"] = ddb.structure.formula
            df["chneut"] = chneut
            for k in ddb_header_keys:
                df[k] = ddb.header[k]

            # Add path to the DDB file.
            if with_path: df["ddb_path"] = ddb.filepath

            df_list.append(df)

        # Concatenate dataframes.
        return dict2namedtuple(df=pd.concat(df_list, ignore_index=True).sort_values(by="site_index"),
                               becs_list=becs_list)

    def anacompare_epsinf(self, ddb_header_keys=None, chneut=1, tol=1e-3, with_path=False, verbose=0):
        r"""
        Compute (eps^\inf) electronic dielectric tensor for all DDBs in the robot and build DataFrame.
        with Voigt indices as columns + metadata. Useful for convergence studies.

        Args:
            ddb_header_keys: List of keywords in the header of the DDB file
                whose value will be added to the Dataframe.
            chneut: Anaddb input variable. See official documentation.
            tol: Elements below this value are set to zero.
            with_path: True to add DDB path to dataframe
            verbose: verbosity level. Set it to a value > 0 to get more information

        Return: ``namedtuple`` with the following attributes::

            df: DataFrame with Voigt indices as columns.
            epsinf_list: List of |DielectricTensor| objects with eps^{inf}
        """
        ddb_header_keys = [] if ddb_header_keys is None else list_strings(ddb_header_keys)
        df_list, epsinf_list = [], []
        for label, ddb in self.items():
            # Invoke anaddb to compute e_inf
            einf, _ = ddb.anaget_epsinf_and_becs(chneut=chneut, verbose=verbose)
            epsinf_list.append(einf)
            df = einf.get_voigt_dataframe(tol=tol)

            # Add metadata to the dataframe.
            df["formula"] = ddb.structure.formula
            df["chneut"] = chneut
            for k in ddb_header_keys:
                df[k] = ddb.header[k]

            # Add path to the DDB file.
            if with_path: df["ddb_path"] = ddb.filepath
            df_list.append(df)

        # Concatenate dataframes.
        return dict2namedtuple(df=pd.concat(df_list, ignore_index=True), epsinf_list=epsinf_list)

    def anacompare_eps0(self, ddb_header_keys=None, asr=2, chneut=1, tol=1e-3, with_path=False, verbose=0):
        """
        Compute (eps^0) dielectric tensor for all DDBs in the robot and build DataFrame.
        with Voigt indices as columns + metadata. Useful for convergence studies.

        Args:
            ddb_header_keys: List of keywords in the header of the DDB file
                whose value will be added to the Dataframe.
            asr, chneut, dipdip: Anaddb input variable. See official documentation.
            tol: Elements below this value are set to zero.
            with_path: True to add DDB path to dataframe
            verbose: verbosity level. Set it to a value > 0 to get more information

        Return: ``namedtuple`` with the following attributes::

            df: DataFrame with Voigt as columns.
            eps0_list: List of |DielectricTensor| objects with eps^0.
            dgen_list: List of DielectricTensorGenerator.
        """
        ddb_header_keys = [] if ddb_header_keys is None else list_strings(ddb_header_keys)
        df_list, eps0_list, dgen_list = [], [], []
        for label, ddb in self.items():
            # Invoke anaddb to compute e_0
            gen = ddb.anaget_dielectric_tensor_generator(asr=asr, chneut=chneut, dipdip=1, verbose=verbose)
            dgen_list.append(gen)
            eps0_list.append(gen.eps0)
            df = gen.eps0.get_voigt_dataframe(tol=tol)

            # Add metadata to the dataframe.
            df["formula"] = ddb.structure.formula
            df["asr"] = asr
            df["chneut"] = chneut
            #df["dipdip"] = dipdip
            for k in ddb_header_keys:
                df[k] = ddb.header[k]

            # Add path to the DDB file.
            if with_path: df["ddb_path"] = ddb.filepath

            df_list.append(df)

        # Concatenate dataframes.
        return dict2namedtuple(df=pd.concat(df_list, ignore_index=True),
                               eps0_list=eps0_list, dgen_list=dgen_list)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        if all(ddb.has_at_least_one_atomic_perturbation() for ddb in self.abifiles):
            print("Invoking anaddb through anaget_phonon_plotters...")
            r = self.anaget_phonon_plotters()
            for fig in r.phbands_plotter.yield_figs(): yield fig
            for fig in r.phdos_plotter.yield_figs(): yield fig

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If ``nbpath`` is None, a temporary file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        anaget_phonon_plotters_kwargs = ( "\n"
            '\tnqsmall=10, ndivsm=20, asr=2, chneut=1, dipdip=1, dos_method="tetra",\n'
            '\tlo_to_splitting=False, ngqpt=None, qptbounds=None,\n'
            '\tanaddb_kwargs=None, verbose=0')

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.DdbRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("""#dfq = robot.get_dataframe_at_qpoint(qpoint=None, units="meV")"""),
            nbv.new_code_cell("r = robot.anaget_phonon_plotters(%s)" % anaget_phonon_plotters_kwargs),
            nbv.new_code_cell("r.phbands_plotter.get_phbands_frame()"),
            nbv.new_code_cell("r.phbands_plotter.ipw_select_plot()"),
            nbv.new_code_cell("r.phdos_plotter.ipw_select_plot()"),
            nbv.new_code_cell("r.phdos_plotter.ipw_harmonic_thermo()"),
        ])

        # Mixins
        nb.cells.extend(self.get_baserobot_code_cells())

        return self._write_nb_nbpath(nb, nbpath)
