# coding: utf-8
"""
Classes and functions for parsing the ONCVPSP output file and plotting the results.
"""
from __future__ import annotations

import os
import re
import tempfile
import numpy as np

from collections import namedtuple
from typing import Union, Any  #, List, Optional
from monty.functools import lazy_property
from monty.collections import AttrDict, dict2namedtuple
from abipy.core.atom import NlkState, RadialFunction, RadialWaveFunction, l2char
from abipy.ppcodes.base_parser import BaseParser


def is_integer(s: Any) -> bool:
    """True if object `s` in an integer."""
    try:
        c = float(s)
        return int(c) == c
    except (ValueError, TypeError):
        return False

# Object returned by self._grep
GrepResults = namedtuple("GrepResults", "data, start, stop")

# Used to store ae and pp quantities (e.g wavefunctions) in a single object.
AePsNamedTuple = namedtuple("AePsNamedTuple", "ae, ps")

ConvData = namedtuple("ConvData", "l energies values")

AtanLogDer = namedtuple("AtanLogDer", "l, energies, values")


class OncvParser(BaseParser):
    """
    Object to read and extract data from the output file of oncvpsp.

    Attributes:
        atsym
        Z
        nc
        nv
        iexc
        psfile

    Example:

        parser = OncvParser(filename)
        parser.scan()

        # To access data:
        parser.radial_wavefunctions

        # To plot data with matplotlib.
        p = parser.get_plotter()
        p.plot_atanlogder_econv()

    """
    # TODO Improve fully-relativistic case.

    def scan(self, verbose: int = 0) -> None:
        """
        Scan the output file, set `run_completed` attribute.

        Raises: self.Error if invalid file.
        """
        try:
            return self._scan(verbose=verbose)
        except Exception as exc:
            raise self.Error(f"Exception while parsing: {self.filepath}") from exc

    def _scan(self, verbose: int = 0) -> None:
        if not os.path.exists(self.filepath):
            raise self.Error(f"File {self.filepath} does not exist")

        # Read data and store it in lines
        self.lines = []
        import io
        with io.open(self.filepath, "rt") as fh:
            for i, line in enumerate(fh):
                #print(line)
                line = line.strip()
                self.lines.append(line)

                if verbose and line.startswith("fcfact*="):
                    print(line)

                if line.startswith("DATA FOR PLOTTING"):
                    self.run_completed = True

                # lines that contain the word ERROR but do not seem to indicate an actual teminating error
                acceptable_error_markers = [
                  'run_config: ERROR for fully non-local  PS atom,'
                ]

                if "ERROR" in line:
                    # Example:
                    # test_data: must have fcfact>0.0 for icmod= 1
                    # ERROR: test_data found   1 errors; stopping
                    if line in acceptable_error_markers:
                        self._warnings.append("\n".join(self.lines[i-1:i+1]))
                    else:
                        self._errors.append("\n".join(self.lines[i-1:i+1]))

                if "WARNING" in line:
                    self._warnings.append("\n".join(self.lines[i:i+2]))

                if "GHOST(+)" in line:
                    self._warnings.append(line)
                if "GHOST(-)" in line:
                    self._errors.append(line)

        # Get gendate, calc_type and version
        # scalar-relativistic version 2.1.1, 03/26/2014
        # scalar-relativistic version 3.0.0 10/10/2014
        toks = self.lines[1].replace(",", " ").split()
        self.gendate = toks.pop(-1)
        self.calc_type, self.version = toks[0], toks[2]

        #print(self.version)
        self.major_version, self.minor_version, self.patch_level = tuple(map(int, self.version.split(".")))[:3]

        # Read configuration (not very robust because we assume the user didn't change the template but oh well)

        # Also, handle atom with z > 100 in which atsym and z are not separated by white space e.g.
        #
        # atsym  z   nc   nv     iexc    psfile
        # Rf104.00   10    8       4      both

        header = "# atsym  z    nc    nv    iexc   psfile"
        for i, line in enumerate(self.lines):

            if line.startswith("# atsym"):
                values = self.lines[i + 1].split()
                if len(values) != 6:
                  # Rf104.00   10    8       4      both
                  l = values.pop(0)
                  atmsym, z = l[0:2], l[2:]
                  values.insert(0, z)
                  values.insert(0, atmsym)
                  #print(values)

                keys = header[1:].split()
                # assert len(keys) == len(values)
                # Store values in self.
                for k, v in zip(keys, values):
                    # Convert nc and nv to int.
                    if k in ("nc", "nv", "iexc"): v = int(v)
                    if k in ("z", ): v = float(v)
                    setattr(self, k, v)
                break

        # Parse pseudization options for the local part.
        header = "# lloc, lpopt,  rc(5),   dvloc0"
        self.rc5 = None
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                tokens = self.lines[i + 1].split()
                #print("tokens", tokens)
                self.lloc = int(tokens[0])
                self.lptopt = int(tokens[1])
                self.rc5 = float(tokens[2])
                self.dvloc0 = float(tokens[3])
                break

        if self.rc5 is None:
            raise self.Error(f"Cannot find magic line starting with `{header}` in: {self.filepath}")

        # Parse ATOM and Reference configuration
        # Example:
        """
        #
        #   n    l    f        energy (Ha)
            1    0    2.00    -6.5631993D+01
            2    0    2.00    -5.1265474D+00
            2    1    6.00    -3.5117357D+00
            3    0    2.00    -3.9736459D-01
            3    1    2.00    -1.4998149D-01

        full rel
        in version 4, there no difference between FR and SR
        in version 3, the FR version has:

        #   n    l    f              energy (Ha)
        #   n    l    f        l+1/2             l-1/2
            1    0    2.00    -2.4703720D+03
            2    0    2.00    -4.2419865D+02

        """
        if self.relativistic and self.major_version <= 3:
            header = "#   n    l    f        l+1/2             l-1/2"
        else:
            header = "#   n    l    f        energy (Ha)"

        nc, nv = self.nc, self.nv

        for i, line in enumerate(self.lines):
            if line.startswith(header):
                beg, core = i + 1, [],
                for c in range(nc):
                    n, l, f = self.lines[beg+c].split()[:3]
                    if is_integer(f):
                        f = str(int(float(f)))
                    else:
                        f = "%.1f" % f
                    core.append(n + l2char[l] + "^%s" % f)
                #self.core = "$" + " ".join(core) + "$"

                beg, valence = i + nc + 1, []
                for v in range(nv):
                    #print("lines[beg+v]", self.lines[beg+v])
                    n, l, f = self.lines[beg+v].split()[:3]
                    if is_integer(f):
                        f = str(int(float(f)))
                    else:
                        #print(f)
                        f = "%.1f" % float(f)

                    valence.append(n + l2char[l] + "^{%s}" % f)

                #self.valence = "$" + " ".join(valence) + "$"
                #print("core", self.core, "valence",self.valence)
                break
        else:
            raise self.Error(f"Cannot find header:\n`{header}`\nin output file {self.filepath}")

    @lazy_property
    def lmax(self) -> int:
        # Read lmax (not very robust because we assume the user didn't change the template but oh well)
        header = "# lmax"
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                return int(self.lines[i+1])
                break
        else:
            raise self.Error(f"Cannot find line with `#lmax` in: {self.filepath}")

    def to_string(self, verbose: int = 0 ) -> str:
        """
        String representation.
        """
        lines = []
        app = lines.append

        if not hasattr(self, "calc_type"):
            app("Object is empty. Call scan method to analyze output file")
            return "\n".join(lines)

        if not self.run_completed:
            app("completed: %s" % self.run_completed)
            return "\n".join(lines)

        app("%s, oncvpsp version: %s, date: %s" % (self.calc_type, self.version, self.gendate))

        from pprint import pformat
        app(pformat(self.get_results()))

        if self.warnings:
            lines.extend(self.warnings)

        if self.errors:
            lines.extend(self.errors)

        return "\n".join(lines)

    def __str__(self) -> str:
        return self.to_string()

    @property
    def relativistic(self) -> bool:
        """True if fully-relativistic calculation."""
        return self.calc_type in ("fully-relativistic", "relativistic")

    @lazy_property
    def rc_l(self) -> dict[int, float]:
        """
        Core radii as a function of l extracted from the output file.
        """
        rc_l = {}
        header = "#   l,   rc,"
        for i, line in enumerate(self.lines):
            if line.startswith(header):
                beg = i + 1
                nxt = 0
                while True:
                    ln = self.lines[beg + nxt]
                    if ln.startswith("#"): break
                    tokens = ln.split()
                    #print("line:", ln, "\ntokens", tokens)
                    l, rc = int(tokens[0]), float(tokens[1])
                    rc_l[l] = rc
                    nxt += 1

        if not rc_l:
            raise self.Error(f"Cannot find magic line starting with `{header}` in: {self.filepath}")

        return rc_l

    @lazy_property
    def kinerr_nlk(self) -> dict[NlkState, namedtuple]:
        """
        Dictionary with the error on the kinetic energy indexed by nlk.
        """

        # In relativist mode we write data inside the following loops:

        #do l1=1,lmax+1
        #   ll=l1-1
        #   if(ll==0) then
        #    mkap=1
        #   else
        #    mkap=2
        #   end if
        #   do ikap=1,mkap
        #       if(ikap==1) then
        #         kap=-(ll+1)
        #       else
        #         kap=  ll
        #       end if

        kinerr_nlk = {}

        if self.major_version > 3:
            # Calculating optimized projector #   1
            #
            #  for l=   0

            re_start = re.compile(r"^Calculating optimized projector #\s+(?P<iproj>\d+)")

        else:
            # Calculating first optimized projector for l=   0
            re_start = re.compile(r"^Calculating (?P<iproj>(first|second)) optimized projector for l=\s+(?P<l>\d+)")
            # TODO: In FR mode, we have
            #Calculating first optimized projector for l=   0
            #Calculating second optimized projector for l=   0

        nlk = None
        iproj_l_seen = set()

        for i, line in enumerate(self.lines):

            m = re_start.match(line)

            if m:
                # Extract iproj and l.
                if self.major_version > 3:
                    # for l=   0
                    iproj = int(m.group("iproj"))
                    l = int(self.lines[i+2].split("=")[-1].strip())
                else:
                    iproj = m.group("iproj")
                    iproj = {"first": 0, "second": 1}[iproj]
                    l = int(m.group("l"))

                k = None
                if self.relativistic:
                    k = 1
                    if (iproj, l) in iproj_l_seen: k= 2
                    iproj_l_seen.add((iproj, l))

                # Use n index to store iprj index.
                nlk = NlkState(n=iproj, l=l, k=k)
                #print("nlk:", nlk)
                continue

            # Now parse the following section associated to nlk

            #Energy error per electron        Cutoff
            #     Ha          eV             Ha
            #     0.01000     0.27211       27.01
            #     0.00100     0.02721       52.82
            #     0.00010     0.00272       66.22
            #     0.00001     0.00027       75.37

            if line.startswith("Energy error per electron        Cutoff"):
                values_ha, ecuts = [], []
                for j in range(4):
                    tokens = self.lines[i+2+j].split()
                    #print("tokens:", tokens)
                    if not tokens: break
                    err_ha, err_ev, ecut = map(float, tokens)
                    values_ha.append(err_ha)
                    ecuts.append(ecut)

                if nlk is None:
                    raise self.Error("Cannot find nlk quantum numbers")

                self._check_nlk_key(nlk, kinerr_nlk, "kinerr_nlk")

                kinerr_nlk[nlk] = dict2namedtuple(ecuts=ecuts, values_ha=values_ha)

        if not kinerr_nlk:
            raise self.Error(f"Cannot parse convergence profile in: {self.filepath}")

        return kinerr_nlk

    @staticmethod
    def _check_nlk_key(nlk, d, dict_name):

        if nlk in d:
            ks = "\n\t".join(str(k) for k in d)
            raise RuntimeError(f"nlk state `{nlk}` is already in {dict_name}:\nKeys:\n\t{ks}")

    @lazy_property
    def potentials(self) -> dict[int, RadialFunction]:
        """
        Dict with radial functions with the non-local and local potentials indexed by l.
        l = -1 corresponds to the local part (if present).
        """
        #radii, charge, pseudopotentials (ll=0, 1, lmax)
        #!p   0.0099448   4.7237412  -7.4449470 -14.6551019
        vl_data = self._grep("!p").data
        lmax = len(vl_data[0]) - 3
        assert lmax == self.lmax

        # From 0 up to lmax
        ionpots_l = {}
        for l in range(lmax + 1):
            ionpots_l[l] = RadialFunction("Ion Pseudopotential, l=%d" % l, vl_data[:, 0], vl_data[:, 2+l])

        # Local part is stored with l == -1 if lloc=4, not present if lloc = l
        vloc = self._grep("!L").data
        if vloc is not None:
            ionpots_l[-1] = RadialFunction("Local part, l=%d" % -1, vloc[:, 0], vloc[:, 1])

        return ionpots_l

    @lazy_property
    def densities(self) -> dict[str, RadialFunction]:
        """
        Dictionary with charge densities on the radial mesh.
        """
        # radii, charge, core charge, model core charge
        # !r   0.0100642   4.7238866  53.4149287   0.0000000
        rho_data = self._grep("!r").data

        return dict(
            rhoV=RadialFunction("Valence charge", rho_data[:, 0], rho_data[:, 1]),
            rhoC=RadialFunction("Core charge", rho_data[:, 0], rho_data[:, 2]),
            rhoM=RadialFunction("Model charge", rho_data[:, 0], rho_data[:, 3])
        )

    @lazy_property
    def radial_wfs(self) -> AePsNamedTuple:
        """
        Read and set the radial wavefunctions for the bound states.

        Usage:

            ae_wfs, ps_wfs = self.radial_wfs.ae, self.radial_wfs.ps

            for nlk, ae_wf in ae_wfs.items():
                ps_wf, l, k = ps_wfs[nlk], nlk.l, nlk.k
        """
        return self._get_radial_wavefunctions(what="bound_states")

    @lazy_property
    def scattering_wfs(self) -> AePsNamedTuple:
        """
        Read and set the scattering wavefunctions.
        """
        return self._get_radial_wavefunctions(what="scattering_states")

    def _get_radial_wavefunctions(self, what: str) -> AePsNamedTuple:
        # For scalar-relativistic bound states, we have
        #
        #    n= 1,  l= 0, all-electron wave function, pseudo w-f
        #
        #    &     0    0.009945   -0.092997    0.015273

        # bound states, fully-relativistic case:
        #
        #    n= 1,  l= 0  kap=-1, all-electron wave function, pseudo w-f
        #
        #    &     0    0.009955    0.066338    0.000979

        # For the scattering states (scalar and relativistic case)
        #
        # scattering, iprj= 2,  l= 1, all-electron wave function, pseudo w-f
        #
        # scattering, iprj= 2,  l= 1, kap= 1, all-electron wave function, pseudo w-f

        ae_waves, ps_waves = {}, {}

        beg = 0
        while True:
            g = self._grep("&", beg=beg)
            if g.data is None: break
            beg = g.stop + 1

            # Get header two lines above.
            header = self.lines[g.start - 2]

            if what == "bound_states":
                if header.startswith("scattering,"):
                    continue
            elif what == "scattering_states":
                if not header.startswith("scattering,"):
                    continue
                header = header.replace("scattering,", "")
            else:
                raise ValueError(f"Invalid value of what: `{what}`")

            #print(header)

            if not self.relativistic:
                # n= 1,  l= 0, all-electron wave function, pseudo w-f
                n, l = header.split(",")[0:2]
                n = int(n.split("=")[1])
                l = int(l.split("=")[1])
                kap = None
            else:
                # n= 1,  l= 0,  kap=-1, all-electron wave function, pseudo w-f
                if self.major_version <= 2: header = header.replace("kap=", ", kap=")
                n, l, kap = header.split(",")[0:3]
                n = int(n.split("=")[1])
                l = int(l.split("=")[1])
                kap = int(kap.split("=")[1])

            nlk = NlkState.from_nlkap(n=n, l=l, kap=kap)
            #print("Got nlk state:", nlk)

            rmesh = g.data[:, 1]
            ae_wf = g.data[:, 2]
            ps_wf = g.data[:, 3]

            self._check_nlk_key(nlk, ae_waves, "ae_waves")

            ae_waves[nlk] = RadialWaveFunction(nlk, str(nlk), rmesh, ae_wf)
            ps_waves[nlk] = RadialWaveFunction(nlk, str(nlk), rmesh, ps_wf)

        return AePsNamedTuple(ae=ae_waves, ps=ps_waves)

    @lazy_property
    def projectors(self) -> dict[NlkState, RadialFunction]:
        """
        Dict with projector wave functions indexed by nlk.
        """
        #
        #@     0    0.009945    0.015274   -0.009284
        beg = 0
        magic = "@"
        if self.major_version > 3: magic = "!J"

        # if(ikap==1) then
        #   write(6,'(a,i6,6(f12.6,1x))') '!J',-ll,rr(ii), &
        #        (vkb(ii,jj,l1,ikap),jj=1,nproj(l1))
        # else
        #   write(6,'(a,i6,6(f12.6,1x))') '!J',ll,rr(ii), &
        #        (vkb(ii,jj,l1,ikap),jj=1,nproj(l1))

        projectors_nlk = {}
        while True:
            g = self._grep(magic, beg=beg)
            if g.data is None: break
            beg = g.stop + 1

            rmesh = g.data[:, 1]
            l = int(g.data[0, 0])

            k = None
            if self.relativistic:
                k = 2
                if l <= 0: k = 1

            for n in range(len(g.data[0]) - 2):
                nlk = NlkState(n=n + 1, l=abs(l), k=k)
                #print("Got projector with: %s" % str(nlk))

                if nlk in projectors_nlk:
                    raise self.Error("nlk state `{nlk}` is already in projectors_nlk")

                projectors_nlk[nlk] = RadialWaveFunction(nlk, str(nlk), rmesh, g.data[:, n + 2])

        return projectors_nlk

    @lazy_property
    def atan_logders(self) -> AePsNamedTuple:
        """
        Atan of the log derivatives for different l-values.
        """
        #log derivativve data for plotting, l= 0
        #atan(r * ((d psi(r)/dr)/psi(r))), r=  1.60
        #l, energy, all-electron, pseudopotential
        #
        #!      0    2.000000    0.706765    0.703758
        ae_atan_logder_l, ps_atan_logder_l = {}, {}

        lstop = self.lmax + 1
        if self.major_version > 3:
            lstop = min(self.lmax + 2, 4)

        if not self.relativistic:
            l_list = list(range(lstop))
        else:
            # Order with (l, -l) for plotting purposes.
            l_list = []
            for l in range(lstop):
                if l != 0:
                    l_list.extend((l, -l))
                else:
                    l_list.append(0)

        for l in l_list:
            tag = "!      %d" % l if l >= 0 else "!     %d" % l
            data = self._grep(tag=tag).data
            if data is None:
                raise self.Error(f"Cannot find logder for l: {l}")
            assert l == int(data[0, 0])

            ae_atan_logder_l[l] = AtanLogDer(l=l, energies=data[:, 1], values=data[:, 2])
            ps_atan_logder_l[l] = AtanLogDer(l=l, energies=data[:, 1], values=data[:, 3])

        return AePsNamedTuple(ae=ae_atan_logder_l, ps=ps_atan_logder_l)

    @lazy_property
    def kene_vs_ecut(self) -> dict[int, ConvData]:
        """
        Dict with the convergence of the kinetic energy versus ecut for different l-values.
        """
        #convergence profiles, (ll=0,lmax)
        #!C     0    5.019345    0.010000
        #...
        #!C     1   19.469226    0.010000
        # TODO: This does not take into account scattering states or n > 1
        conv_l = {}

        for l in range(self.lmax + 1):
            data = self._grep(tag="!C     %d" % l).data
            conv_l[l] = ConvData(l=l, energies=data[:, 1], values=data[:, 2])

        return conv_l

    @lazy_property
    def hints(self) -> dict:
        """
        Hints for the cutoff energy as provided by oncvpsp.
        """
        # Extract the hints
        hints = 3 * [-np.inf]
        for i in range(3):
            for l in range(self.lmax + 1):
                hints[i] = max(hints[i], self.kene_vs_ecut[l].energies[-i-1])
        hints.reverse()

        # Truncate to the nearest int
        hints = [np.rint(h) for h in hints]
        # print("hints:", hints)

        return dict(
            low={"ecut": hints[0], "pawecutdg": hints[0]},
            normal={"ecut": hints[1], "pawecutdg": hints[1]},
            high={"ecut": hints[2], "pawecutdg": hints[2]}
        )

    def get_results(self) -> AttrDict:
        """
        Return the most important results extracted from the output file.
        """
        # Init return values
        #d = AttrDict(
        #    max_ecut=None,
        #    max_atan_logder_l1err=None,
        #    max_psexc_abserr=None,
        #    herm_err=None,
        #    nwarns=len(self.warnings)
        #    nerrs=len(self.errors)
        #)

        # Get the max ecut estimated by oncvpsp.
        # TODO: Should take into account scattering states.
        max_ecut = max(self.kene_vs_ecut[l].energies[-1] for l in self.kene_vs_ecut)

        # Compute the l1 error in atag(logder) between AE and PS
        from scipy.integrate import cumtrapz
        max_l1err = 0.0
        for l in self.atan_logders.ae:
            f1, f2 = self.atan_logders.ae[l], self.atan_logders.ps[l]

            abs_diff = np.abs(f1.values - f2.values)
            integ = cumtrapz(abs_diff, x=f1.energies) / (f1.energies[-1] - f1.energies[0])
            max_l1err = max(max_l1err, integ[-1])

        # Read Hermiticity error and compute the max value of PSP excitation error=
        # Hermiticity error    4.8392D-05
        # PSP excitation error=  1.56D-10
        herm_tag, pspexc_tag = "Hermiticity error", "PSP excitation error="
        herm_err, max_psexc_abserr = None, -np.inf

        for line in self.lines:
            i = line.find(herm_tag)
            if i != -1:
                herm_err = float(line.split()[-1].replace("D", "E"))

            i = line.find(pspexc_tag)
            if i != -1:
                max_psexc_abserr = max(max_psexc_abserr, abs(float(line.split()[-1].replace("D", "E"))))

        return AttrDict(
            max_ecut=max_ecut,
            max_atan_logder_l1err=max_l1err,
            max_psexc_abserr=max_psexc_abserr,
            herm_err=herm_err,
            nwarns=len(self.warnings),
            nerrs=len(self.errors),
        )

    def find_string(self, s: str) -> int:
        """
        Returns the index of the first line containing string s.
        Raises self.Error if s cannot be found.
        """
        for i, line in enumerate(self.lines):
            if s in line:
                return i
        else:
            raise self.Error(f"Cannot find `{s}` in lines")

    def get_input_str(self) -> str:
        """String with the ONCVPSP input file."""
        try:
            # oncvpsp 3.2.3
            i = self.find_string("<INPUT>")
            j = self.find_string("</INPUT>")
            return "\n".join(self.lines[i+1:j]) + "\n"
        except self.Error:
            # oncvpsp => 4
            i = self.find_string("Reference configufation results")
            return "\n".join(self.lines[:i])

    def get_psp8_str(self) -> Union[str, None]:
        """
        Return string with the pseudopotential data in psp8 format.
        Return None if field is not present.
        """
        start, stop = None, None
        for i, line in enumerate(self.lines):
            if 'Begin PSPCODE8' in line: start = i
            if start is not None and 'END_PSP' in line:
                stop = i
                break

        if start is None and stop is None: return None
        ps_data = "\n".join(self.lines[start+1:stop])

        if "<INPUT>" not in ps_data:
            # oncvpsp <= 3.2.2 --> Append the input to ps_data (note XML markers)
            # oncvpsp >= 3.2.3 --> Input is already there
            ps_data += "\n\n<INPUT>\n" + self.get_input_str() + "</INPUT>\n"

        return ps_data

    def get_upf_str(self) -> Union[str, None]:
        """
        Return string with the pseudopotential data in upf format.
        Return None if field is not present.
        """
        start, stop = None, None
        for i, line in enumerate(self.lines):
            if "Begin PSP_UPF" in line: start = i
            if start is not None and 'END_PSP' in line:
                stop = i
                break

        if start is None and stop is None: return None
        return "\n".join(self.lines[start+1:stop])

    def get_plotter(self): # -> Union[OncvPlotter, None]:
        """
        Return an instance of OncvPlotter or None
        """
        from abipy.ppcodes.oncv_plotter import OncvPlotter
        try:
            return OncvPlotter(self)
        except Exception as exc:
            print(exc)
            #raise
            return None

    def _grep(self, tag: str, beg: int = 0) -> GrepResults:
        """
        Finds the first field in the file with the specified tag.
        `beg` gives the initial position in the file.
        """
        data, stop, intag = [], None, -1

        if beg >= len(self.lines):
            raise ValueError(f"beg {beg} > len(lines) ({len(self.lines)})")

        for i, l in enumerate(self.lines[beg:]):
            l = l.lstrip()
            if l.startswith(tag):
                if intag == -1:
                    intag = beg + i
                data.append([float(c) for c in l.split()[1:]])
            else:
                # Exit because we know there's only one section starting with 'tag'
                if intag != -1:
                    stop = beg + i
                    break

        if not data:
            return GrepResults(data=None, start=intag, stop=stop)
        else:
            return GrepResults(data=np.array(data), start=intag, stop=stop)

    def gnuplot(self) -> None:
        """
        Plot the results with gnuplot.
        Based on the `replot.sh` script provided by the oncvpsp code.
        """
        outfile = self.filepath
        base = os.path.basename(outfile)
        gnufile = base + ".scr"
        plotfile = base + ".plot"
        temp = base + ".tmp"

        from monty.os import cd
        from subprocess import check_call
        workdir = tempfile.mkdtemp()
        print(f"Working in: {workdir}")

        with cd(workdir):
            check_call("awk 'BEGIN{out=0};/GNUSCRIPT/{out=0}; {if(out == 1) {print}}; \
                                /DATA FOR PLOTTING/{out=1}' %s > %s" % (outfile, plotfile), shell=True)

            check_call("awk 'BEGIN{out=0};/END_GNU/{out=0}; {if(out == 1) {print}}; \
                                /GNUSCRIPT/{out=1}' %s > %s" % (outfile, temp), shell=True)

            check_call('sed -e 1,1000s/t1/"%s"/ %s > %s' % (plotfile, temp, gnufile), shell=True)

            try:
                check_call(["gnuplot", gnufile])
            except KeyboardInterrupt:
                print("Received KeyboardInterrupt")

        os.rmdir(workdir)
