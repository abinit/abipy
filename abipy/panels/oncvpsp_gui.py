#!/usr/bin/env python
from __future__ import annotations

import sys
import os
import functools
import itertools
import param
import panel as pn
import pandas as pd
import panel.widgets as pnw

from monty.termcolor import cprint
from monty.string import list_strings
from abipy.panels.core import AbipyParameterized, depends_on_btn_click, mpl, dfc, ButtonContext, Loading
from abipy.ppcodes.ppgen import OncvGenerator
from abipy.ppcodes.oncv_parser import OncvParser


GE_ANNOTATED = """
# Copied from ~oncvpsp/doc/32_Ge_annotated.dat

# ATOM AND REFERENCE CONFIGURATION
#
### atsym  atomic symbol
### z  atomic number
### nc number of core states
### nv number of valence states
### iexc  exchange-correlation functional: 1-Wigner, 2-Hedin-Lundquist,
###  3-Perdew-Zunger-Ceperly-Alder, 4-Perdew-Burke-Enzerhof
### psfile format of pseudopotential file, psp8 for ABINIT, upf for PWSCF
#
# atsym, z, nc, nv, iexc   psfile
    Ge  32.0   5   3   3   psp8
#
#
### n principal quantum number
### l angular momentum
### f occupancy (MUST be >0)
#
# n, l, f  (nc+nv lines)
    1    0    2.0
    2    0    2.0
    2    1    6.0
    3    0    2.0
    3    1    6.0
    3    2   10.0
    4    0    2.0
    4    1    2.0
#
# PSEUDOPOTENTIAL AND OPTIMIZATION
#
###lmax maximum angular momentum for which psp is calculated (<=3)
#
# lmax
    2
#
#
### l angular momentum
### rc  core radius for this l
### ep  energy at which psp is generated (eigenvalue inserted for occupied
###     state in reference configuration, positive energy must be specified
###     for barrier-confined "scattering" state for unoccupied l <=lmax
###     A small positive energy is usually  good (0.1-0.25 Ha).
### ncon number of constraints for pseudo wave function to match all-
###      electron wave function at rc, value + ncon-1 derivatives,
###      must be between 3 and 5 ("M" in the paper, Eq.(6))
### nbas number of basis functions.  Must be between ncon+2 and ncon+5
###      ("N" in the paper, Eqs.(4-5))
### qcut wave vector defining "residual energy" in the RRKJ method
###      ("q_c" in the paper, Eq.(1)
#
# l, rc, ep, ncon, nbas, qcut  (lmax+1 lines, l's must be in order)
    0    2.60   -0.00    4    8    5.00
    1    2.60   -0.00    4    8    5.20
    2    2.00    0.00    4    9    8.40
#
# LOCAL POTENTIAL
#
### lloc angular momentum whose semi-local psp is taken as local.  lloc=4
###      denotes a smooth polynomial continuation of the all-electron
###      potential.  If lloc<=lmax, remaining data are ignored, but
###      must be there (zeros are OK).  The rc corresponding to lloc
###      MUST BE THE MINIMUM rc for lloc<=lmax, or <= the minumum for
###      lloc=4 (usually the best choice)
### lpopt type of polynomial continuation for lloc=4. values 1-5
###       permitted.
###     1) match 2 derivatives, r^0,2,4
###     2) match 2 derivatives, r^0,4,6
###     3) match 3 derivatives, r^0,4,5,6
###     4) match 3 derivatives, r^0,4,6,8
###     5) match 3 derivatives, r^0,2,4,6
### dvloc0 shift of r=0 potential from basic continuation (Ha) for lloc=4
###        depends on lpopt above as follows, with x=(r/rc)
###        1) dvloc0*(1-x^2)^3
###        2) dvloc0*(1-x^4)^3
###        3-5) dvloc0*(1-x^4)^4
#
# lloc, lpopt, rc(5), dvloc0
    4    5    2.0    0.0
#
# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs
#
### l angular momentum
### nproj number of projectors, 1 or 2.  automatically set to 0 for l=lloc
### debl energy  added to basic psp  energy ep for 2nd projector
###      automatically reset to match 2nd projector with 2nd bound state
###      at this l when it is occupied (ie., the psp is generated for a
###      correcponding-l shallow core state)
#
# l, nproj, debl  (lmax+1 lines, l's in order)
    0    2    1.50
    1    2    1.50
    2    2    1.50
#
# MODEL CORE CHARGE
#
### icmod 0: no non-linear core correction charge. 1: smooth monotonic
###       polynomial model core charge fit at "matching" rc following
###       reference 35. For icmod = 2, 3, 4 see doc/core_correction.txt
### fcfact radius for above determined by  rho_core(r)=fcfact*rho_pseudo_
###        valence(r) values 0.25-0.5 are usually good (look at plots)
###       For icmod = 3, fcfact has a different meaning and a third
###       artument rcfact is added (see core_correction.txt)
###       For icmod = 4,, fcfact is ignored but must be present (0.0 OK)
#
# icmod, fcfact, (rcfact)
    0    0.25
#
# LOG DERIVATIVE ANALYSIS
#
### epsh1 lower energy limit for "phase-shift-like" log-derivative plots
###       should be below the nearest core level for this l to make sure
###       there are no ghosts, but -2.0 Ha usually is OK  When semi-cores
###       are treated as valence, this should be below the lowest core
###       energy
### epsh2 upper energy limit, 2.0 usually good
### depsh energy mesh interval for plot, 0.02 usually good enough
#
# epsh1, epsh2, depsh
   -2.0  2.0  0.02
#
# OUTPUT GRID
#
### rlmax maximum radius for Abinit psp code 8 format output.  must be
###       greater than maximum rc (including lloc=4 rc), but also determines
###       range of diagnostic plots, so  ~2-3*rcmax is usually good
### drl mesh spacing of linear radial mesh for Abinit output. 0.02 is good
###     for "softer" psps, 0.01 is probably better with 1st row, 3d's, or
###     semi-core psps.
#
# rlmax, drl
    4.0  0.01
#
# TEST CONFIGURATIONS
#
### ncnf number of test configurations (<=4)  The reference config is always
###      run first as a consistency check.  core always is the reference core
###      configuration.  The excitation energy of the configuration for the
###      pseudo-atom will be compared with the all-electron result.
# ncnf
    4
#
### n principal quantum number for all-electron atom
### l angular momentum
### f occupancy
#
#   nvcnf (repeated ncnf times)
#         number of valence states in this test configuration
#   n, l, f  (nvcnf lines, repeated follwing nvcnf's ncnf times)
#   n    l    f
    3
    3    2   10.00
    4    0    1.00
    4    1    2.00
#
    3
    3    2   10.00
    4    0    2.00
    4    1    1.00
#
    3
    3    2   10.00
    4    0    1.00
    4    1    1.00
#
    3
    3    2   10.00
    4    0    1.00
    4    2    1.00
#
"""

class Lparams(AbipyParameterized):
    """
    Stores all the oncvpsp pseudization parameters for a given l.
    """

    l = param.Integer(None, bounds=(0, None))
    rc = param.Number(None, bounds=(0, None))
    ep = param.Number(None, bounds=(None, None))
    ncon = param.Integer(None, bounds=(0, 9))
    nbas = param.Integer(None, bounds=(0, 9))
    qcut = param.Number(None, bounds=(0, None))
    nproj = param.Integer(None, bounds=(1, 5,))
    debl = param.Number(None, bounds=(None, None))



class Nlf(AbipyParameterized):
    """
    Stores the value of n, l and occupancy f.
    """
    n = param.Integer(None, bounds=(1, None))
    l = param.Integer(None, bounds=(0, None))
    f = param.Number(None, bounds=(0.0, None))



class OncvInput(AbipyParameterized):
    """
    A parametrized class with all the Oncvpsp input variables,
    usually constructed from an external file.
    """

    atsym = param.String()
    z = param.Integer(1, bounds=(1, None))
    nc = param.Integer(0, bounds=(0, None))
    nv = param.Integer(0, bounds=(0, None))
    iexc = param.Integer(0)
    psfile = param.String("both")

    nlf_list = param.List() #[], item_type=Nlf)
    lmax = param.Integer(0, bounds=(0, 4))
    lparams = param.List() #[], item_type=Lparams)

    lloc = param.Integer(4)
    lpopt = param.Integer(5)
    rc5 = param.Number(0.0)
    dvloc0 = param.Number(0.0)

    icmod = param.Integer(0, bounds=(0, 3))
    fcfact = param.Number(0.0, bounds=(0, None))
    rcfact = param.Number(0.0, bounds=(0, None))

    epsh1 = param.Number(-25)
    epsh2 = param.Number(+25)
    depsh = param.Number(0.02)
    rlmax = param.Number(6.0)
    drl = param.Number(0.1)

    nconf = param.Integer(0, bounds=(0, None))

    @classmethod
    def from_file(cls, path: str) -> OncvInput:
        """
        Initialize the object from file.
        """
        with open(path, "rt") as fh:
            return cls.from_string(fh.read())

    @classmethod
    def from_string(cls, string: str) -> OncvInput:
        """Build the object from a string."""

        def parse_line(line: str, types: str):
            """
            Helper function to parse line, convert each token according to
            character in `types` and return list of values or value dependening on num char in types.
            s -> str, i -> int, f -> float.
            """
            try:
                i = line.find("#")
                if i != -1:
                    line = line[:i]
                tokens = line.split()

                if len(tokens) != len(types):
                    raise ValueError(f"Expecting {len(types)} tokens but got {len(tokens)} tokens: `{tokens}`")

                outs = []
                for tok, typ in zip(tokens, types):
                    if typ == "s":
                        outs.append(str(tok))
                    elif typ == "i":
                        fval = float(tok)
                        ival = int(fval)
                        if fval != ival:
                            raise TypeError(f"Expecting int in line {line}, got : `{tok}`")
                        outs.append(ival)
                    elif typ == "f":
                        outs.append(float(tok))
                    else:
                        raise TypeError(f"Invalid type: `{typ}`")

                return outs if len(outs) > 1 else outs[0]

            except Exception as exc:
                print(exc)
                raise ValueError(f"Invalid line: `{line}`, expecting types: `{types}`")

        lines = [l for l in string.split("\n") if not l.startswith("#")]
        lines = [l for l in lines if l.strip()]

        # atsym z nc nv iexc psfile
        atsym, z, nc, nv, iexc, psfile = parse_line(lines.pop(0), "siiiis")
        psfile = "both"

        nlf_list = []
        for i in range(nc + nv):
            # n l f
            n, l, f = parse_line(lines.pop(0), "iif")
            nlf_list.append(Nlf(n=n, l=l, f=f))

        lmax = parse_line(lines.pop(0), "i")

        params_l = {}
        for i in range(lmax + 1):
            # l rc ep ncon nbas qcut
            l, rc, ep, ncon, nbas, qcut = parse_line(lines.pop(0), "iffiif")
            assert i == l
            params_l[l] = dict(l=l, rc=rc, ep=ep, ncon=ncon, nbas=nbas, qcut=qcut)

        # lloc lpopt rc5 dvloc0
        lloc, lpopt, rc5, dvloc0 = parse_line(lines.pop(0), "iiff")

        for l in range(lmax + 1):
            # l nproj debl
            l, nproj, debl = parse_line(lines.pop(0), "iif")
            params_l[l].update(nproj=nproj, debl=debl)

        lparams = [Lparams(**d) for d in params_l.values()]

        # icmod fcfact [rcfact]
        l = lines.pop(0)
        rcfact = 0.0
        try:
            icmod, fcfact, rcfact = parse_line(l, "iff")
        except ValueError:
            icmod, fcfact = parse_line(l, "if")

        # epsh1 epsh2 depsh
        epsh1, epsh2, depsh = parse_line(lines.pop(0), "fff")
        # rlmax drl
        rlmax, drl = parse_line(lines.pop(0), "ff")
        #nconf = parse_line(lines.pop(0), "i")
        nconf = 0

        locs = locals()
        d = {k: locs[k]for k in [
            "atsym",
            "z",
            "nc",
            "nv",
            "iexc",
            "psfile",
            "nlf_list",
            "lparams",
            "lmax",
            "lloc",
            "lpopt",
            "rc5",
            "dvloc0",
            "icmod",
            "fcfact",
            "rcfact",
            "epsh1",
            "epsh2",
            "depsh",
            "rlmax",
            "drl",
            "nconf",
        ]}

        return cls(**d)


    def __init__(self, **params):
        super().__init__(**params)

    def __str__(self) -> str:
        """String with the ONCVPSP input file."""
        lines = []
        app = lines.append

        app("# ATOM AND REFERENCE CONFIGURATION")
        app("# atsym z nc nv iexc psfile")
        app(f"{self.atsym} {self.z} {self.nc} {self.nv} {self.iexc} {self.psfile}")
        app("# n l f (nc+nv lines)")
        for i, nlf in enumerate(self.nlf_list):
            tag = "# end core" if i + 1 == self.nc else ""
            app(f"{nlf.n} {nlf.l} {nlf.f} {tag}")

        app("# PSEUDOPOTENTIAL AND OPTIMIZATION")
        app("# lmax")
        app(f"{self.lmax}")
        app("# l rc ep ncon nbas qcut (lmax+1 lines, l's must be in order)")
        for p in self.lparams:
            app(f"{p.l} {p.rc:.2f} {p.ep:.2f} {p.ncon} {p.nbas} {p.qcut:.2f}")
        app("# LOCAL POTENTIAL")
        app("# lloc lpopt rc5 dvloc0")
        app(f"{self.lloc} {self.lpopt} {self.rc5:.2f} {self.dvloc0}")
        app("# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs")
        app("# l nproj debl (lmax+1 lines, l's in order")
        for p in self.lparams:
            app(f"{p.l} {p.nproj} {p.debl:.2f}")

        app("# MODEL CORE CHARGE")
        app("# icmod fcfact rcfact")
        app(f"{self.icmod} {self.fcfact:.2f} {self.rcfact:.2f}")
        app("# LOG DERIVATIVE ANALYSIS")
        app("# epsh1 epsh2 depsh")
        app(f"{self.epsh1} {self.epsh2} {self.depsh}")
        app("# OUTPUT GRID")
        app("# rlmax drl")
        app(f"{self.rlmax} {self.drl}")
        app("# TEST CONFIGURATIONS")
        app("# ncnf")
        app(f"{self.nconf}")

        return "\n".join(lines)

    def get_min_rc(self) -> float:
        """Return the minimum of rc(l) over l."""
        return min(p.rc for p in self.lparams)

    def find_lparam(self, l: int, what: str) -> Tuple(int, float):
        for i, p in enumerate(self.lparams):
            if p.l == l:
                return i, getattr(p, what)

        raise ValueError(f"Cannot find l: {l} in input with what: {what}")


def run_psgen(psgen: OncvGenerator, data: dict) -> dict:
    """
    Helper function to run oncvpsp, parse the outputfiles and return dictionary with results

    Args:
        psgen (OncvGenerator): [description]
        data (dict): input dictionary with extra data to be added to the output dict
    """
    # Init return values
    max_ecut = None
    max_atan_logder_l1err = None
    max_psexc_abserr = None
    herm_err = None
    status = None
    nwarns = None
    nerrs = None

    try:
        psgen.start()
        retcode = psgen.wait()

        results = psgen.parser.get_results()

        if results is not None:
            max_ecut = results.max_ecut
            max_atan_logder_l1err = results.max_atan_logder_l1err
            max_psexc_abserr = results.max_psexc_abserr
            herm_err = results.herm_err
            status = str(psgen.status)
            nwarns = len(psgen.parser.warnings)
            nerrs = len(psgen.parser.errors)

    except Exception as exc:
        print("Exception in run_psgen", exc)

    d =  dict(
        status=status,
        max_ecut=max_ecut,
        max_atan_logder_l1err=max_atan_logder_l1err,
        max_psexc_abserr=max_psexc_abserr,
        herm_err=herm_err,
        nwarns=nwarns,
        nerrs=nerrs,
    )

    data = data.copy()
    data.update(**d)
    return data



def build_mesh(x0: float, num: int, step: float, direction: str) -> list:
    """
    Generate a linear mesh of step `step` that is centered on x0 if
    directions == "centered" or a mesh that starts/ends at x0 if direction is `>`/`<`.
    """

    if direction == "centered":
        start = x0 - num * step
        return [start + i * step for i in range(2 * num + 1)]
    elif direction in (">", "<"):
        start = x0
        if direction == "<": step = -abs(step)
        return sorted([start + i * step for i in range(num)])
    else:
        raise ValueError(f"Invalid direction: `{direction}`")



class OncvGui(AbipyParameterized):

    calc_type = param.ObjectSelector(default="scalar-relativistic",
                                     objects=["scalar-relativistic", "fully-relativistic", "non-relativistic"],
                                     label="Relativistic effects")

    max_nprocs =  param.Integer(max(os.cpu_count() // 2, 1), bounds=(1, None))

    dpi = param.Integer(82, bounds=(24, None))

    #in_filepath = param.String("", doc="The path to the oncvps input file.")

    qcut_num =  param.Integer(2, bounds=(1, None))
    qcut_step = param.Number(1, bounds=(0, None))
    qcut_dir = param.Selector(["centered", ">", "<"])

    rc_num =  param.Integer(2, bounds=(1, None))
    rc_step = param.Number(0.1, bounds=(0, None))
    rc_dir = param.Selector(["centered", ">", "<"])

    debl_num =  param.Integer(2, bounds=(1, None))
    debl_step = param.Number(1.0, bounds=(0, None))
    debl_dir = param.Selector(["centered", ">", "<"])

    rc5_num =  param.Integer(2, bounds=(1, None))
    rc5_step = param.Number(0.1, bounds=(0, None))
    rc5_dir = param.Selector(["<", "centered", ">"])

    dvloc0_num =  param.Integer(2, bounds=(1, None))
    dvloc0_step = param.Number(0.5, bounds=(0, None))
    dvloc0_dir = param.Selector(["centered", "<", ">"])

    fcfact_num =  param.Integer(1, bounds=(1, None))
    fcfact_step = param.Number(0.2, bounds=(0, None))
    fcfact_dir = param.Selector(["centered", ">", "<"])

    rcfact_num =  param.Integer(1, bounds=(1, None))
    rcfact_step = param.Number(0.2, bounds=(0, None))
    rcfact_dir = param.Selector(["centered", ">", "<"])

    ace_theme = param.ObjectSelector(default="chrome",
                                     objects=pnw.Ace.param.theme.objects,
                                     doc="Theme of the editor")

    history_idx = param.Integer(default=-1, label="History index")

    @classmethod
    def from_file(cls, path: str) -> OncvGui:
        """
        Build an instance from a file with the oncvpsp input variables.
        """
        return cls(oncv_input=OncvInput.from_file(path), in_filepath=path)

    def __init__(self, oncv_input, in_filepath="", **params):
        super().__init__(**params)

        self.ace_kwargs = dict(sizing_mode='stretch_both', print_margin=False, language='text', height=600,
                          theme="chrome",
                          #theme="dracula",
                          #max_length=150,
                          )

        self.input_ace = pnw.Ace(value=str(oncv_input), **self.ace_kwargs)

        # Add annotated example for documentation purposes.
        self.annotated_example = pn.pane.HTML(f"<pre><code> {GE_ANNOTATED} </code></pre>")

        self.in_filepath = pnw.TextInput(value=in_filepath, placeholder='Enter the filepath...')

        self.out_area = pn.Column("## Oncvpsp output:",  sizing_mode="stretch_width")
        #self.out_runtests = pn.Column("## Basic tests:", sizing_mode="stretch_width")

        # Define buttons
        self.execute_btn = pnw.Button(name="Execute", button_type='primary')
        self.execute_btn.on_click(self.on_execute_btn)

        # This is the directory used to run oncvpsp when the user clicks execute_btn
        #self._execute_stdout_path = None

        # List storing all the inputs.
        self.input_history = []
        self.history_btn = pnw.Button(name="Compare", button_type='primary')
        #self.history_btn.on_click(self.on_history_btn)

        self.rc_qcut_btn = pnw.Button(name="Execute", button_type='primary')

    @param.depends("ace_theme")
    def change_ace_theme(self):
        print("Changing theme")
        self.input_ace.theme = self.ace_theme

    def get_oncv_input(self) -> OncvInput:
        """
        Take the string from the ACE editor and build an oncv input
        with the last changes done by the user
        """
        return OncvInput.from_string(self.input_ace.value)

    def set_oncv_input_string(self, new_string: str) -> OncvInput:
        """
        Update the string in the ACE editor, store previous string in history.
        Return: OncvInput instance.
        """
        #print("Updating input with new_string:\n", new_string)
        self.input_history.append(self.input_ace.value)
        self.input_ace.value = new_string

        return self.get_oncv_input()

    def starmap(self, func, list_of_args):
        import time
        time_start = time.time()

        # Don't use more procs than tasks.
        _max_nprocs_ = self.max_nprocs
        _max_nprocs_ = min(_max_nprocs_, len(list_of_args))

        if _max_nprocs_ == 1:
            values = [func(*args) for args in list_of_args]
        else:
            # TODO: Optimize, use better algo
            # This Pool uses threads instead of multiprocessing
            # Cannot use multiprocessing because we have side effects in psgen
            from multiprocessing.dummy import Pool
            with Pool(processes=self.max_nprocs) as pool:
                values = pool.starmap(func, list_of_args)

        print(f"Done {len(list_of_args)} tasks in {time.time() - time_start:.2f} [s] with {_max_nprocs_} processe(s)")
        return values

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets"""
        d = {}

        main = pn.Column(
            pn.Row(
                self.pws_col(["calc_type", "max_nprocs",
                              "dpi", "ace_theme", "execute_btn"]),
                self.input_ace,
            ),
            pn.Card(self.annotated_example, title='Annotated example', collapsed=True,
                    header_color="blue", sizing_mode="stretch_width"),
            sizing_mode="stretch_width",
        )

        d["Main"] =  pn.Column(main,
                               #pn.layout.Divider(),
                               self.out_area,
                               #self.out_runtests,
                               sizing_mode="stretch_width",
                               )

        d["Rc_qcut_opt"] = self.get_rc_qcut_opt_view()
        d["History"] = self.get_history_view()

        if as_dict: return d
        template = self.get_template_from_tabs(d, template=kwargs.get("template", None))
        #self.tabs = template
        #print(self.tabs)
        return template

    def get_history_view(self):
        return pn.Row(
            self.pws_col(["## History",
                          "history_idx",
                          "history_btn",
                         ]),
            self.on_history_btn
        )

    @depends_on_btn_click('history_btn')
    def on_history_btn(self):
        print("hello")
        hist_len = len(self.input_history)
        idx = self.history_idx
        if hist_len == 0 or idx >= hist_len:
            return pn.Column(f"hist_len == {hist_len}")

        ace_hist = pnw.Ace(value=self.input_history[idx], **self.ace_kwargs)

        fromlines = self.input_history[idx].splitlines()
        tolines = self.input_ace.value.splitlines()
        from difflib import HtmlDiff
        html_table = HtmlDiff().make_table(fromlines, tolines) #, fromdesc='', todesc='', context=False, numlines=5)

        return pn.Column(
                f"## Input at history index: {idx}",
                ace_hist,
                pn.pane.HTML(html_table),
        )

    def get_rc_widgets(self, oncv_input):
        """Return widgets to change the value of rc(l)"""
        menu_items = [(f"l = {l}", str(l)) for l in range(oncv_input.lmax + 1)]
        menu_button = pnw.MenuButton(name='Change rc(l)', items=menu_items, button_type='primary')
        menu_button.on_click(self.on_change_rc)
        rc_l = {p.l: p.rc for p in oncv_input.lparams}
        help_str = f"""
Here one can change the value of rc(l).

The present values of rc_l are: {rc_l}
"""
        return pn.WidgetBox(menu_button,
                            *[self.param[k] for k in ("rc_num", "rc_step", "rc_dir")],
                            help_str)

    def get_qcut_widgets(self, oncv_input):
        """Return widgets to change the value of qc(l)"""
        menu_items = [(f"l = {l}", str(l)) for l in range(oncv_input.lmax + 1)]
        menu_button = pnw.MenuButton(name='Change qcut(l)', items=menu_items, button_type='primary')
        menu_button.on_click(self.on_change_qcut)
        qc_l = {p.l: p.qcut for p in oncv_input.lparams}
        help_str = f"""
Here one can change the value of qcut(l).

The present values are: {qc_l}
"""
        return pn.WidgetBox(menu_button,
                            *[self.param[k] for k in ("qcut_num", "qcut_step", "qcut_dir")],
                            help_str)

    def get_debl_widgets(self, oncv_input):
        """Return widgets to change the value of debl(l)"""
        menu_items = [(f"l = {l}", str(l)) for l in range(oncv_input.lmax + 1)]
        menu_button = pnw.MenuButton(name='Change debl(l)', items=menu_items, button_type='primary')
        menu_button.on_click(self.on_change_debl)
        help_str = f"""
Here one can change the value of debl(l) with fixed nproj(l).
"""
        return pn.WidgetBox(menu_button,
                            *[self.param[k] for k in ("debl_num", "debl_step", "debl_dir")],
                            help_str)

    def get_rc5_widgets(self, oncv_input):
        """Return widgets to change the value of rc5"""
        btn = pnw.Button(name="Run", button_type='primary')
        btn.on_click(self.on_change_rc5)
        help_str = f"""
Here one can change the value of rc5 for vloc.

The present value of rc5 is {oncv_input.rc5} and min(rc) is: {oncv_input.get_min_rc()}
"""
        return pn.WidgetBox(*[self.param[k] for k in ("rc5_num", "rc5_step", "rc5_dir")],
                            btn,
                            help_str)

    def get_dvloc0_widgets(self, oncv_input):
        """Return widgets to change the value of dvloc0"""
        btn = pnw.Button(name="Run", button_type='primary')
        btn.on_click(self.on_change_dvloc0)
        help_str = f"""
Here one can change the value of dvloc0 for vloc.

The present value of dvloc0 is {oncv_input.dvloc0} with lpopt: {oncv_input.lpopt}
"""
        return pn.WidgetBox(*[self.param[k] for k in ("dvloc0_num", "dvloc0_step", "dvloc0_dir")],
                            btn,
                            help_str)

    def get_rhomodel_widgets(self, oncv_input):
        """Return widgets to change the parameters for the model core charge"""
        btn = pnw.Button(name="Run", button_type='primary')
        btn.on_click(self.on_change_rhomodel)
        help_str = f"""
Here one can change the parameters for the model core charge.

The present value of icmod is {oncv_input.icmod} with fcfact: {oncv_input.fcfact} and rcfact: {oncv_input.rcfact}
"""
        keys = ("rcfact_num", "rcfact_step", "rcfact_dir",
                "fcfact_num", "fcfact_step", "fcfact_dir")
        return pn.WidgetBox(*[self.param[k] for k in keys],
                            btn,
                            help_str)

    def gridplot_psgens(self, psgens, titles, func_names="plot_atanlogder_econv"):
        """
        Return a GridBox with the figures obtained by calling `plotter.func_name`
        for all the PseudoGenerators in psgens.
        """
        # Generate plots with titles by calling `func_name`.
        _m = functools.partial(mpl, with_divider=False, dpi=self.dpi)
        func_names = list_strings(func_names)
        figs = []
        for psgen, title in zip(psgens, titles):
            plotter = psgen.parser.get_plotter()
            if plotter is not None:
                for func_name in func_names:
                    plot_func = getattr(plotter, func_name)
                    figs.append(_m(plot_func(show=False, title=title, fig_close=True)))

        # Insert the figures in a GridBox.
        nfigs = len(figs)
        nrows, ncols = 1, 1
        if nfigs > 1:
            ncols = 2 if len(func_names) == 1 else len(func_names)
            nrows = nfigs // ncols + nfigs % ncols
        elif nfigs == 0:
            nrows, ncols = 0, 0

        return pn.GridBox(*figs, ncols=ncols, nrows=nrows)

    def on_change_qcut(self, event):
        """
        Change the value of qcut(l), run oncvpsp and show the results.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            # Get initial qc(l) from input.
            l = int(event.new)
            oncv_input = self.get_oncv_input()
            i0, qcut0 = oncv_input.find_lparam(l, "qcut")

            # Define list of qc to be tested and build list of OncvGenerator.
            qcut_values = build_mesh(qcut0, self.qcut_num, self.qcut_step, self.qcut_dir)
            psgens = []

            try:
                for qcut in qcut_values:
                    oncv_input.lparams[i0].qcut = qcut
                    psgens.append(OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type))
            finally:
                # Restore the initial value.
                oncv_input.lparams[i0].qcut = qcut0

            # Run each generator, get results, build pandas Dataframe and show it in out_area.
            tasks = [(psgen, {"qcut": qcut}) for psgen, qcut in zip(psgens, qcut_values)]
            d_list = self.starmap(run_psgen, tasks)

            dfw = self._build_table(tasks, d_list)

            head = pn.Row(pn.Column(
                            f"## Qcut optimization for l: {l}. Click the icon to update the input",
                            dfw
                            ),
                         self.get_qcut_widgets(oncv_input))

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add plots
            grid = self.gridplot_psgens(psgens, [f"qc = {qc:.2f}" for qc in qcut_values])
            col.append(grid)

            self.out_area.objects = col.objects

    def _build_table(self, tasks, d_list):
        """
        Build and return a Tabulator widget with the results.
        Also, register callbacks so that it is possible to
        update the input file by clicking on the icon in the last column.
        """

        df = pd.DataFrame(d_list, columns=list(d_list[0].keys()))

        # Sort results by max_ecut and add buttons to trigger callbacks
        df = df.sort_values("max_ecut")
        dfw = pn.widgets.Tabulator(df, buttons={'accept': '<i class="fa fa-print"></i>'})

        def update_input(event):
            #print(f'Clicked {event.column!r} on row {event.row}')
            # Use the index to get the psgen as we have sorted along max_ecut
            idx = dfw.value.index[event.row]
            #print("select idx:", idx)
            psgen = tasks[idx][0]
            if event.column == "accept" and psgen.status == psgen.S_OK:
                oncv_input = self.set_oncv_input_string(psgen.input_str)
                self._update_out_area(psgen, oncv_input)

        dfw.on_click(update_input)
        return dfw

    def on_change_debl(self, event):
        """
        Change the value of debl(l), run oncvpsp and show the results.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            # Get initial qc(l) from input.
            l = int(event.new)
            oncv_input = self.get_oncv_input()

            # Define list of qc to be tested and build list of OncvGenerator.
            i0, debl0 = oncv_input.find_lparam(l, "debl")
            debl_values = build_mesh(debl0, self.debl_num, self.debl_step, self.debl_dir)
            psgens = []

            try:
                for debl in debl_values:
                    oncv_input.lparams[i0].debl = debl
                    psgens.append(OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type))
            finally:
                # Restore the initial value.
                oncv_input.lparams[i0].debl = debl0

            # Run each generator, get results, build pandas Dataframe and show it in out_area.
            tasks = [(psgen, {"debl": debl}) for psgen, debl in zip(psgens, debl_values)]
            d_list = self.starmap(run_psgen, tasks)

            dfw = self._build_table(tasks, d_list)

            head = pn.Row(pn.Column(
                            f"## Debl optimization for l: {l}. Click the icon to update the input",
                            dfw,
                            ),
                         self.get_debl_widgets(oncv_input))

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add plots
            grid = self.gridplot_psgens(psgens, [f"debl = {debl:.2f}" for debl in debl_values],
                                        func_names=["plot_atan_logders"])
            col.append(grid)

            self.out_area.objects = col.objects

    def on_change_rc5(self, event):
        """
        Change the value of rc5 for the local part, run oncvpsp and show the results.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            oncv_input = self.get_oncv_input()

            # Define list of qc to be tested and build list of OncvGenerator.
            rc5 = oncv_input.rc5
            rc5_values = build_mesh(rc5, self.rc5_num, self.rc5_step, self.rc5_dir)
            psgens = []

            try:
                for new_rc in rc5_values:
                    oncv_input.rc5 = new_rc
                    psgens.append(OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type))
            finally:
                # Restore the initial value.
                oncv_input.rc5 = rc5

            # Run each generator, get results, build pandas Dataframe and show it in out_area.
            tasks = [(psgen, {"rc5": rc5}) for psgen, rc5 in zip(psgens, rc5_values)]
            d_list = self.starmap(run_psgen, tasks)

            dfw = self._build_table(tasks, d_list)

            head = pn.Row(pn.Column(
                            f"## Rc5 optimization. Click the icon to update the input",
                            dfw,
                            ),
                         self.get_rc5_widgets(oncv_input))

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add plots
            grid = self.gridplot_psgens(psgens, [f"rc5 = {rc5:.2f}" for rc5 in rc5_values],
                                        func_names=["plot_atanlogder_econv", "plot_potentials"])
            col.append(grid)

            self.out_area.objects = col.objects

    def on_change_dvloc0(self, event):
        """
        Change the value of dvloc0 for the local part, run oncvpsp and show the results.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            oncv_input = self.get_oncv_input()

            # Define list of qc to be tested and build list of OncvGenerator.
            dvloc0 = oncv_input.dvloc0
            dvloc_values = build_mesh(dvloc0, self.dvloc0_num, self.dvloc0_step, self.dvloc0_dir)
            psgens = []

            try:
                for new_dvloc in dvloc_values:
                    oncv_input.dvloc0 = new_dvloc
                    psgens.append(OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type))
            finally:
                # Restore the initial value.
                oncv_input.dvloc0 = dvloc0

            # Run each generator, get results, build pandas Dataframe and show it in out_area.
            tasks = [(psgen, {"dvloc0": dvloc}) for psgen, dvloc in zip(psgens, dvloc_values)]
            d_list = self.starmap(run_psgen, tasks)

            dfw = self._build_table(tasks, d_list)

            head = pn.Row(pn.Column(
                            f"## dvloc0 optimization. Click the icon to update the input",
                            dfw,
                            ),
                         self.get_dvloc0_widgets(oncv_input))

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add plots
            grid = self.gridplot_psgens(psgens, [f"dvloc = {dvloc:.2f}" for dvloc in dvloc_values],
                                        func_names=["plot_atanlogder_econv", "plot_potentials"])
            col.append(grid)

            self.out_area.objects = col.objects

    def on_change_rc(self, event):
        """
        Change the value of rc(l), run oncvpsp and show the results.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            # Get initial rc(l) from input.
            l = int(event.new)
            oncv_input = self.get_oncv_input()

            # Define list of qc to be tested and build list of OncvGenerator.
            i0, rc0 = oncv_input.find_lparam(l, "rc")
            rc_values = build_mesh(rc0, self.rc_num, self.rc_step, self.rc_dir)
            psgens = []

            try:
                for rc in rc_values:
                    oncv_input.lparams[i0].rc = rc
                    psgens.append(OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type))
            finally:
                # Restore the initial value.
                oncv_input.lparams[i0].rc = rc0

            # Run each generator, get results, build pandas Dataframe and show it in out_area.
            tasks = [(psgen, {"rc": rc}) for psgen, rc in zip(psgens, rc_values)]
            d_list = self.starmap(run_psgen, tasks)

            dfw = self._build_table(tasks, d_list)

            head = pn.Row(pn.Column(
                            f"## Rc optimization for l: {l}. Click the icon to update the input",
                            dfw,
                            ),
                         self.get_rc_widgets(oncv_input),
                         )

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add plots
            grid = self.gridplot_psgens(psgens, [f"rc = {rc:.2f} for l: {l}" for rc in rc_values])
            col.append(grid)

            self.out_area.objects = col.objects

    def on_change_rhomodel(self, event):
        """
        Change the parameters for the model core charge, run oncvpsp and show the results.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            # Get initial values from input.
            oncv_input = self.get_oncv_input()
            icmod = oncv_input.icmod
            if icmod == 0: return
            fcfact0 = oncv_input.fcfact
            rcfact0 =  oncv_input.rcfact

            # Define list of values to be tested and build list of OncvGenerator.
            fcfact_values = build_mesh(fcfact0, self.fcfact_num, self.fcfact_step, self.fcfact_dir)
            rcfact_values = build_mesh(rcfact0, self.rcfact_num, self.rcfact_step, self.rcfact_dir)

            psgens = []
            tasks = []
            titles = []

            if icmod in (3, 4):
                # fcfact x rcfact
                for fc, rc in itertools.product(fcfact_values, rcfact_values):
                    try:
                        oncv_input.fcfact = fc
                        oncv_input.rcfact = rc
                        psgens.append(OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type))
                        tasks.append((psgens[-1], {"fcfact": fc, "rcfact": rc}))
                        titles.append(f"fcfact: {fc}, rcfact: {rc}")
                    finally:
                        # Restore the initial value.
                        oncv_input.fcfact = fcfact0
                        oncv_input.rcfact = rcfact0

            elif icmod in (1, 2):
                # Only fcfact is used here.
                for fc in fcfact_values:
                    try:
                        oncv_input.fcfact = fc
                        psgens.append(OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type))
                        tasks.append((psgens[-1], {"fcfact": fc}))
                        titles.append(f"fcfact: {fc}")
                    finally:
                        # Restore the initial value.
                        oncv_input.fcfact = fcfact0
            else:
                raise ValueError(f"Invalid icmod: {icmod}")

            # Run each generator, get results, build pandas Dataframe and show it in out_area.
            d_list = self.starmap(run_psgen, tasks)

            dfw = self._build_table(tasks, d_list)

            head = pn.Row(pn.Column(
                            f"## Rho model optimization for icmod: {icmod}. Click the icon to update the input",
                            dfw,
                            ),
                         self.get_rhomodel_widgets(oncv_input),
                         )

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add plots
            grid = self.gridplot_psgens(psgens, titles,
                                       func_names=["plot_densities", "plot_den_formfact"])
            col.append(grid)

            self.out_area.objects = col.objects

    def get_rc_qcut_opt_view(self):
        oncv_input = self.get_oncv_input()

        menu_items = [(f"l = {l}", str(l)) for l in range(oncv_input.lmax + 1)]
        menu_button = pnw.MenuButton(name='Change qcut(l)', items=menu_items, button_type='primary')
        menu_button.on_click(self.on_change_rc_qcut)
        qc_l = {p.l: p.qcut for p in oncv_input.lparams}
        rc_l = {p.l: p.rc for p in oncv_input.lparams}
        help_str = f"""
Here one can change the value of rc(l) and qcut(l).

The present values of qc_l are: {qc_l}
The present values of rc_l are: {rc_l}
"""
        self.rc_qcut_out_area = pn.Column(sizing_mode="stretch_width")
        wbox = pn.WidgetBox(menu_button,
                            *[self.param[k] for k in ("qcut_num", "qcut_step", "qcut_dir")],
                            *[self.param[k] for k in ("rc_num", "rc_step", "rc_dir")],
                             help_str)

        return pn.Row(wbox, self.rc_qcut_out_area, sizing_mode="stretch_width")

    #@depends_on_btn_click('rc_qcut_btn')
    def on_change_rc_qcut(self, event):
        """
        Generate pseudos using a grid of (rc, qcut) values for given l.
        """
        with ButtonContext(event.obj), Loading(self.rc_qcut_out_area):
            # Get initial rc(l) from input.
            l = int(event.new)
            oncv_input = self.get_oncv_input()

            # Define list of qc to be tested and build list of OncvGenerator.
            i0, rc0 = oncv_input.find_lparam(l, "rc")
            rc_values = build_mesh(rc0, self.rc_num, self.rc_step, self.rc_dir)

            # Define list of qc to be tested and build list of OncvGenerator.
            i0, qcut0 = oncv_input.find_lparam(l, "qcut")
            qcut_values = build_mesh(qcut0, self.qcut_num, self.qcut_step, self.qcut_dir)

            def rq_prod():
                return itertools.product(rc_values, qcut_values)

            psgens = []
            try:
                for rc, qcut in rq_prod():
                    oncv_input.lparams[i0].rc = rc
                    oncv_input.lparams[i0].qcut = qcut
                    psgens.append(OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type))
            finally:
                # Restore the initial value.
                oncv_input.lparams[i0].rc = rc0
                oncv_input.lparams[i0].qcut = qcut0

            # Run each generator, get results, build pandas Dataframe and show it in out_area.
            tasks = [(psgen, {"rc": rc, "qcut": qcut}) for psgen, (rc, qcut) in zip(psgens, rq_prod())]
            d_list = self.starmap(run_psgen, tasks)

            dfw = self._build_table(tasks, d_list)

            head = pn.Row(pn.Column(
                            f"## Rc/qcut optimization for l: {l}. Click the icon to update the input",
                            dfw,
                            ),
                         )

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add plots:
            grid = self.gridplot_psgens(psgens, [f"rc = {rc:.2f}, qc = {qc:.2f} for l: {l}"
                                        for rc, qc in rq_prod()])
            col.append(grid)

            self.rc_qcut_out_area.objects = col.objects

    #@depends_on_btn_click('execute_btn')
    def on_execute_btn(self, event):
        """
        Build a new generator from the input file, run it and update out_area.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            oncv_input = self.get_oncv_input()
            psgen = OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type)
            print("Running in workdir:", psgen.workdir)
            psgen.start()
            retcode = psgen.wait()
            #psget.check_status()

            #if psgen.status != psgen.S_OK:
            #    cprint("oncvpsp returned %s. Exiting" % retcode, "red")
            #    return 1

            #out_path = self._execute_stdout_path = psgen.stdout_path

            # Parse the output file
            #onc_parser = OncvParser(out_path)
            #onc_parser.scan()
            #if not onc_parser.run_completed:
            #    cprint("oncvpsp output is not complete. Exiting", "red")
            #    return 1

            ## Build the plotter and add figures to out_area.
            #plotter = onc_parser.get_plotter()

            # TODO:
            # Tranfer final output file.

            self._update_out_area(psgen, oncv_input)

    def _update_out_area(self, psgen, oncv_input):

        with Loading(self.out_area):
            #self.psgen_to_save = psgen
            plotter = psgen.parser.get_plotter()

            if plotter is None:
                self.out_area.objects = pn.Column("## Plotter is None")
                return

            _m = functools.partial(mpl, with_divider=False, dpi=self.dpi)

            save_btn = pnw.Button(name="Save output", button_type='primary')
            save_btn.on_click(self.on_save_btn)

            new_rows = [
                pn.Row(save_btn, self.in_filepath),
                pn.layout.Divider(),
                "## Pseudized Wavefunctions",
                pn.Row(_m(plotter.plot_radial_wfs(show=False)),
                       self.get_rc_widgets(oncv_input)),
                pn.Row(_m(plotter.plot_radial_wfs(what="scattering_states", show=False))),
                pn.layout.Divider(),
                "## Logder and convergence profile",
                pn.Row(_m(plotter.plot_atanlogder_econv(show=False)),
                       self.get_qcut_widgets(oncv_input),
                       self.get_debl_widgets(oncv_input)),
                pn.layout.Divider(),
                "## Pseudized local part",
                pn.Row(_m(plotter.plot_potentials(show=False)),
                       self.get_rc5_widgets(oncv_input),
                       self.get_dvloc0_widgets(oncv_input)),
                pn.layout.Divider(),
                "## Model core charge",
                pn.Row(_m(plotter.plot_densities(show=False)),
                       self.get_rhomodel_widgets(oncv_input)),
                pn.Row(_m(plotter.plot_den_formfact(show=False))),
                pn.layout.Divider(),
                "## Projectors",
                pn.Row(_m(plotter.plot_projectors(show=False))),
                pn.layout.Divider(),
            ]

            self.out_area.objects = new_rows
            #self.tabs[0].active = 1

    def on_save_btn(self, event):
        with ButtonContext(event.obj):
            print("on_save_button")
            #self._execute_stdout_path
            out_path = self.in_filepath.value
            if not out_path:
                raise ValueError("out_path cannot be empty.")

            #self.psgen_to_save = psgen

            # TODO: Transfer final output file.
            #shutil.copy(psgen.stdout_path, out_path)
