#!/usr/bin/env python
from __future__ import annotations

import sys
import os
import functools
import param
import panel as pn
import pandas as pd
import panel.widgets as pnw

from monty.termcolor import cprint
from abipy.panels.core import AbipyParameterized, depends_on_btn_click, mpl, dfc, ButtonContext, Loading
from pseudo_dojo.ppcodes.ppgen import OncvGenerator
from pseudo_dojo.ppcodes.oncvpsp import OncvOutputParser, MultiPseudoGenDataPlotter
#from pseudo_dojo.core.dojoreport import DojoReport
#from pseudo_dojo.refdata.nist import database as nist


#    # atsym z nc nv iexc psfile
#    Si 14 1 4 4 psp8
#    # n l f
#    1 0 2.0
#    2 0 2.0
#    2 1 6.0
#    3 0 2.0
#    3 1 2.0
#    # lmax
#    2
#    # l rc ep ncon nbas qcut
#    0 1.6 0.0 3 8 9.0
#    1 0.90 0.0 3 8 14.0
#    2 1.9 0.05 3 8 5.0
#    # lloc lpopt rc5 dvloc0
#    4 5 0.9 0.0
#    # l nproj debl
#    0 2 2.5
#    1 2 2.5
#    2 2 1.5
#    icmod fcfact rcfact
#    3 4.0 1.3
#    # epsh1 epsh2 depsh
#    -12.0 12.0 0.02
#    # rlmax drl
#    6.0 0.01
#     0

class Lparams(AbipyParameterized):
    """Stores all the oncvpsp pseudization parameters for a given l."""

    l = param.Integer(None, bounds=(0, None), doc="")
    rc = param.Number(None, bounds=(0, None), doc="")
    ep = param.Number(None, bounds=(0, None), doc="")
    ncon = param.Integer(None, bounds=(0, 9), doc="") # ??
    nbas = param.Integer(None, bounds=(0, 9), doc="") # ??
    qcut = param.Number(None, bounds=(0, None), doc="")
    nproj = param.Integer(None, bounds=(1, 5,), doc="")
    debl = param.Number(None, bounds=(None, None), doc="")

    #def __panel__(self):
    #    return pn.Column(
    #            f"# l: {self.l}",
    #            self.param["rc"],
    #            self.param["ep"],
    #            self.param["ncon"],
    #            self.param["nbas"],
    #            self.param["qcut"],
    #            self.param["nproj"],
    #            self.param["debl"],
    #        )


class Nlf(AbipyParameterized):
    """
    Stores the value of n, l and occupancy f.
    """

    n = param.Integer(None, bounds=(1, None), doc="")
    l = param.Integer(None, bounds=(0, None), doc="")
    f = param.Number(None, bounds=(0.0, None), doc="")

    #def __panel__(self):
    #    s = f"## n = {self.n}, l = {self.l}"
    #    return self.pws_row([s, "f"])


class OncvInput(AbipyParameterized):
    """
    A parametrized class with all the Oncvpsp input variables,
    typicall constructed from an external file.
    """

    atsym = param.String(doc="")
    z = param.Integer(1, bounds=(1, None), doc="")
    nc = param.Integer(0, bounds=(0, None), doc="")
    nv = param.Integer(0, bounds=(0, None), doc="")
    iexc = param.Integer(0, doc="")
    psfile = param.String("both", doc="")

    nlf_list = param.List() #[], item_type=Nlf)
    lmax = param.Integer(0, bounds=(0, 4), doc="")
    lparams = param.List() #[], item_type=Lparams)

    lloc = param.Integer(4, doc="")
    lpopt = param.Integer(5, doc="")
    rc5 = param.Number(0.0, doc="")
    dvloc0 = param.Number(0.0, doc="")

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
                tokens = line.rstrip("#").split()
                if len(tokens) != len(types):
                    raise ValueError(f"Expecting {len(types)} tokens while got {len(tokens)} tokens: {tokens}")

                outs = []
                for tok, typ in zip(tokens, types):
                    if typ == "s":
                        outs.append(str(tok))
                    elif typ == "i":
                        outs.append(int(tok))
                    elif typ == "f":
                        outs.append(float(tok))
                    else:
                        raise TypeError(f"Invalid type: `{typ}`")

                return outs if len(outs) > 1 else outs[0]

            except Exception as exc:
                print(exc)
                raise ValueError(f"Invalid line: `{line}`, expecting types: `{types}`")

        lines = [l for l in string.split("\n") if not l.startswith("#")]

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
        app("# asym z nc nv iexc psfile")
        app(f"{self.atsym} {self.z} {self.nc} {self.nv} {self.iexc} {self.psfile}")
        app("# n l f")
        for nlf in self.nlf_list:
            app(f"{nlf.n} {nlf.l} {nlf.f}")
        app("# lmax")
        app(f"{self.lmax}")
        app("# l rc ep ncon nbas qcut")
        for p in self.lparams:
            app(f"{p.l} {p.rc} {p.ep} {p.ncon} {p.nbas} {p.qcut}")
        app("# lloc lpopt rc5 dvloc0")
        app(f"{self.lloc} {self.lpopt} {self.rc5} {self.dvloc0}")
        app("# l nproj debl")
        for p in self.lparams:
            app(f"{p.l} {p.nproj} {p.debl}")
        app("# icmod fcfact rcfact")
        app(f"{self.icmod} {self.fcfact} {self.rcfact}")
        app("# epsh1 epsh2 depsh")
        app(f"{self.epsh1} {self.epsh2} {self.depsh}")
        app("# rlmax drl")
        app(f"{self.rlmax} {self.drl}")
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

    try:
        psgen.start()
        retcode = psgen.wait()

        if psgen.results is None and psgen.status == psgen.S_OK:
            psgen.check_status()

        if psgen.results is not None:
            max_ecut = psgen.results.max_ecut
            max_atan_logder_l1err = psgen.results.max_atan_logder_l1err
            max_psexc_abserr = psgen.results.max_psexc_abserr
            herm_err = psgen.results.herm_err
            status=str(psgen.status),
            nwarns = len(psgen.warnings)

    except Exception as exc:
        print("Exception in run_psgen", exc)

    d =  dict(
        status=status,
        max_ecut=max_ecut,
        max_atan_logder_l1err=max_atan_logder_l1err,
        max_psexc_abserr=max_psexc_abserr,
        herm_err=herm_err,
        nwarns=nwarns,
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



class OncvPanel(AbipyParameterized):

    calc_type = param.ObjectSelector(default="scalar-relativistic",
                                     objects=["scalar-relativistic", "fully-relativistic", "non-relativistic"],
                                     doc="Relativistic effects")

    #max_nprocs =  param.Integer(2, bounds=(1, None), doc="")
    max_nprocs =  param.Integer(max(os.cpu_count() // 2, 1), bounds=(1, None), doc="")

    qcut_num =  param.Integer(2, bounds=(1, None), doc="")
    qcut_step = param.Number(0.5, bounds=(0, None), doc="")
    qcut_dir = param.Selector(["centered", ">", "<"])

    rc_num =  param.Integer(2, bounds=(1, None), doc="")
    rc_step = param.Number(0.1, bounds=(0, None), doc="")
    rc_dir = param.Selector(["centered", ">", "<"])

    debl_num =  param.Integer(2, bounds=(1, None), doc="")
    debl_step = param.Number(1.0, bounds=(0, None), doc="")
    debl_dir = param.Selector(["centered", ">", "<"])

    rc5_num =  param.Integer(2, bounds=(1, None), doc="")
    rc5_step = param.Number(0.1, bounds=(0, None), doc="")
    rc5_dir = param.Selector(["<", "centered", ">"])

    @classmethod
    def from_file(cls, path: str) -> OncvPanel:
        return cls(oncv_input=OncvInput.from_file(path))

    def __init__(self, oncv_input, **params):
        super().__init__(**params)

        self.editor = pn.widgets.Ace(value=str(oncv_input), sizing_mode='stretch_both')
                                     #print_margin=False, language='shell', max_length=150) # height=300,

        self.out_area = pn.Column(sizing_mode="stretch_both")

        # Define buttons
        self.execute_btn = pnw.Button(name="Execute", button_type='primary')
        self.execute_btn.on_click(self.on_execute_btn)

        self.accept_btn = pnw.Button(name="Accept", button_type='primary')
        self.accept_btn.on_click(self.on_accept_btn)

    def get_oncv_input(self) -> OncvInput:
        """
        Take the string from the ACE editor and build an oncv input
        with the last changes done by the user
        """
        return OncvInput.from_string(self.editor.value)

    def starmap(self, func, list_of_args):
        import time
        time_start = time.time()
        # Don't use more procs than tasks.
        max_nprocs_ = self.max_nprocs
        max_nprocs_ = min(max_nprocs_, len(list_of_args))
        if max_nprocs_ == 1:
            values = [func(*args) for args in list_of_args]
        else:
            # TODO: Optimize, use better algo
            # This Pool uses threads instead of multiprocessing
            # Cannot use multiprocessing because we have side effects in psgen"
            from multiprocessing.dummy import Pool
            with Pool(processes=self.max_nprocs) as pool:
                values = pool.starmap(func, list_of_args)

        print(f"Done {len(list_of_args)} tasks in {time.time() - time_start} [s] with {max_nprocs_} processe(s)")
        return values

    def __panel__(self):
        head = pn.Row(
            self.pws_col(["calc_type", "max_nprocs", "execute_btn", "accept_btn"]),
            self.editor,
        )
        return pn.Column(head,
                         #pn.layout.Divider(),
                         self.out_area,
                         sizing_mode="stretch_both"
                         )

    def get_rc_widgets(self, oncv_input):
        """Return widgets used to change the value of rc(l)"""
        menu_items = [(f"l = {l}", str(l)) for l in range(oncv_input.lmax + 1)]
        menu_button = pnw.MenuButton(name='Change rc(l)', items=menu_items, button_type='primary')
        menu_button.on_click(self.on_change_rc)
        help_str = f"""
Here one can change the value of rc(l).

The present value of rc5 is {oncv_input.rc5} and min(rc) is: {oncv_input.get_min_rc()}
"""
        return pn.WidgetBox(menu_button,
                            *[self.param[k] for k in ("rc_num", "rc_step", "rc_dir")],
                            help_str)

    def get_qcut_widgets(self, oncv_input):
        """Return widgets used to change the value of qc(l)"""
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
        """Return widgets used to change the value of debl(l)"""
        menu_items = [(f"l = {l}", str(l)) for l in range(oncv_input.lmax + 1)]
        menu_button = pnw.MenuButton(name='Change debl(l)', items=menu_items, button_type='primary')
        menu_button.on_click(self.on_change_debl)
        help_str = f"""
Here one can change the value of debl(l) with fixed nproj(k).
"""
        return pn.WidgetBox(menu_button,
                            *[self.param[k] for k in ("debl_num", "debl_step", "debl_dir")],
                            help_str)


    def get_rc5_widgets(self, oncv_input):
        """Return widgets used to change the value of rc5"""
        btn = pnw.Button(name="Run", button_type='primary')
        btn.on_click(self.on_change_rc5)
        help_str = f"""
Here one can change the value of rc5 for vloc.

The present value of rc5 is {oncv_input.rc5} and min(rc) is: {oncv_input.get_min_rc()}
"""
        return pn.WidgetBox(*[self.param[k] for k in ("rc5_num", "rc5_step", "rc5_dir")],
                            btn,
                            help_str)

    def gridplot_psgens(self, psgens, titles, func_name="plot_atanlogder_econv"):
        """
        Return a GridBox with the figures obtained by calling `plotter.func_name`
        for all the PseudoGenerators in psgens.
        """
        # Generate plots with tiles by calling `func_name`.
        _m = functools.partial(mpl, with_divider=False)
        figs = []
        for psgen, title in zip(psgens, titles):
            plotter = psgen.plotter
            if plotter is not None:
                plot_func = getattr(psgen.plotter, func_name)
                figs.append(_m(plot_func(show=False, title=title)))

        # Insert the figures in a GridBox.
        nfigs = len(figs)
        nrows, ncols = 1, 1
        if nfigs > 1:
            ncols = 2
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

            df = pd.DataFrame(d_list, columns=list(d_list[0].keys()))
            df.reset_index(drop=True, inplace=True)

            head = pn.Row(pn.Column(
                            f"## Qcut optimization for l: {l}",
                            dfc(df, with_export_btn=False, with_divider=False),
                            ),
                         self.get_qcut_widgets(oncv_input))

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add logder plots:
            grid = self.gridplot_psgens(psgens, [f"qc = {qc:.2f}" for qc in qcut_values])
            col.append(grid)

            self.out_area.objects = col.objects

    def on_change_debl(self, event):
        """
        Change the value of debl(l), run oncvpsp and show the results.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            # Get initial qc(l) from input.
            l = int(event.new)
            oncv_input = self.get_oncv_input()
            i0, debl0 = oncv_input.find_lparam(l, "debl")

            # Define list of qc to be tested and build list of OncvGenerator.
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
            tasks = [(psgen, {"debl": debl}) for psgen, delb in zip(psgens, debl_values)]
            d_list = self.starmap(run_psgen, tasks)

            df = pd.DataFrame(d_list, columns=list(d_list[0].keys()))
            df.reset_index(drop=True, inplace=True)

            head = pn.Row(pn.Column(
                            f"## Debl optimization for l: {l}",
                            dfc(df, with_export_btn=False, with_divider=False),
                            ),
                         self.get_debl_widgets(oncv_input))

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add logder plots:
            grid = self.gridplot_psgens(psgens, [f"debl = {debl:.2f}" for debl in debl_values],
                                        func_name="plot_atan_logders")
            col.append(grid)

            self.out_area.objects = col.objects

    def on_change_rc5(self, event):
        """
        Change the value of rc5 for local part, run oncvpsp and show the results.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            oncv_input = self.get_oncv_input()
            rc5 = oncv_input.rc5

            # Define list of qc to be tested and build list of OncvGenerator.
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

            df = pd.DataFrame(d_list, columns=list(d_list[0].keys()))
            df.reset_index(drop=True, inplace=True)

            head = pn.Row(pn.Column(
                            f"## Rc5 optimization",
                            dfc(df, with_export_btn=False, with_divider=False),
                            ),
                         self.get_rc5_widgets(oncv_input))

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add logder plots:
            grid = self.gridplot_psgens(psgens, [f"rc5 = {rc5:.2f}" for rc5 in rc5_values])
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
            i0, rc0 = oncv_input.find_lparam(l, "rc")

            # Define list of qc to be tested and build list of OncvGenerator.
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
            rows = self.starmap(run_psgen, tasks)

            df = pd.DataFrame(rows, columns=list(rows[0].keys()))
            df.reset_index(drop=True, inplace=True)

            head = pn.Row(pn.Column(
                            f"## Rc optimization for l: {l}",
                            dfc(df, with_export_btn=False, with_divider=False),
                            ),
                         self.get_rc_widgets(oncv_input),
                         )

            col = pn.Column(head, sizing_mode="stretch_width")

            # Add logder plots:
            grid = self.gridplot_psgens(psgens, [f"rc = {rc:.2f} for l: {l}" for rc in rc_values])
            col.append(grid)

            self.out_area.objects = col.objects

    #@depends_on_btn_click('execute_btn')
    def on_execute_btn(self, event):
        """
        Build a new generator from the input file, and add it to the queue.
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            oncv_input = self.get_oncv_input()
            self.out_area.objects = []
            psgen = OncvGenerator(input_str=str(oncv_input), calc_type=self.calc_type)
            print("Running in workdir:", psgen.workdir)
            psgen.start()
            retcode = psgen.wait()
            #psget.check_status()

            #if psgen.status != psgen.S_OK:
            #    cprint("oncvpsp returned %s. Exiting" % retcode, "red")
            #    return 1

            out_path = psgen.stdout_path

            # Parse the output file
            #onc_parser = OncvOutputParser(out_path)
            #onc_parser.scan()
            #if not onc_parser.run_completed:
            #    cprint("oncvpsp output is not complete. Exiting", "red")
            #    return 1

            ## Build the plotter and add figures to out_area.
            #plotter = onc_parser.make_plotter()

            plotter = psgen.plotter
            _m = functools.partial(mpl, with_divider=False)

            rows = [
                pn.Row(_m(plotter.plot_radial_wfs(show=False)),
                       self.get_rc_widgets(oncv_input)),
                pn.layout.Divider(),
                pn.Row(_m(plotter.plot_atanlogder_econv(show=False)),
                       self.get_qcut_widgets(oncv_input),
                       self.get_debl_widgets(oncv_input)),
                pn.layout.Divider(),
                pn.Row(_m(plotter.plot_potentials(show=False)),
                       self.get_rc5_widgets(oncv_input)),
                pn.layout.Divider(),
                pn.Row(_m(plotter.plot_densities(show=False)),
                       _m(plotter.plot_den_formfact(show=False))),
                pn.layout.Divider(),
                pn.Row(_m(plotter.plot_projectors(show=False))),
                pn.layout.Divider(),
            ]

            self.out_area.objects = pn.Column(*rows).objects

    def on_accept_btn(self, event):
        """
        Accept the input file, rerun oncvps with the last input and ...
        """
        with ButtonContext(event.obj), Loading(self.out_area):
            print("In on_accept_btn")
            oncv_input = self.get_oncv_input()
            self.out_area.objects = []
            #self.out_area.objects = pn.Column(*rows).objects


def main():
    #if options.seaborn:
    # Use seaborn settings.
    import seaborn as sns
    sns.set(context="paper", style='darkgrid', palette='deep',
            font='sans-serif', font_scale=1, color_codes=False, rc=None)

    pn.extension('ace')

    app = OncvPanel.from_file(sys.argv[1])
    Template = pn.template.FastListTemplate
    #pn.template.FastListTemplate
    #pn.template.BootstrapTemplate
    #pn.template.GoldenLayoutTemplate
    #pn.template.ReactTemplate(main
    app = Template(main=app.__panel__(), title="Oncvpsp GUI")

    serve_kwargs = dict(debug=True)
    return pn.serve(app, **serve_kwargs)


if __name__ == "__main__":
    main()
