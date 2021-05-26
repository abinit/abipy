""""GUIs for structure."""

import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.panels.core import HasStructureParams, dfc, mpl, ply, depends_on_btn_click


class StructurePanel(HasStructureParams):
    """
    Panel with widgets to interact with an AbiPy Structure
    """

    def __init__(self, structure,  **params):
        self._structure = structure

        # Convert widgets.
        self.output_format = pnw.Select(name="format", value="abinit",
                                        options="abinit,cif,xsf,poscar,qe,siesta,wannier90,cssr,json".split(","))

        # Spglib widgets
        self.spglib_symprec = pnw.Spinner(name="symprec", value=0.01, start=0.0, end=None, step=0.01)
        self.spglib_angtol = pnw.Spinner(name="angtol", value=5, start=0.0, end=None, step=1)

        # K-path widgets
        self.kpath_format = pnw.Select(name="format", value="abinit", options=["abinit", "siesta", "wannier90"])
        self.line_density = pnw.Spinner(name="line density", value=10, step=5, start=0, end=None)
        self.plot_kpath = pnw.Checkbox(name='Plot k-path', value=False)

        # MP-match
        self.mp_match_btn = pnw.Button(name="Connect to Materials Project", button_type='primary')

        # MP-search
        #mp_search_btn = pnw.Button(name="Connect to Materials Project", button_type='primary')
        #mp_api_key

        # GS input generator widgets.
        self.gs_input_btn = pnw.Button(name="Generate input", button_type='primary')
        self.gs_type = pnw.Select(name="GS type", value="scf", options=["scf", "relax"])
        self.kppra = pnw.Spinner(name="kppra", value=1000, step=1000, start=0, end=None)

        self.label2mode = {
            "unpolarized": 'unpolarized',
            "polarized": 'polarized',
            "anti-ferromagnetic": "afm",
            "non-collinear with magnetism": "spinor",
            "non-collinear, no magnetism": "spinor_nomag",
        }

        self.spin_mode = pnw.Select(name="SpinMode", value="unpolarized", options=list(self.label2mode.keys()))

        super().__init__(**params)

    @property
    def structure(self):
        return self._structure

    @pn.depends("output_format.value")
    def convert(self):
        """Convert the input structure to one of the format selected by the user."""
        s = self.structure.convert(fmt=self.output_format.value)
        return self.html_with_clipboard_btn(f"<pre> {s} </pre>")

    @pn.depends("spglib_symprec.value", "spglib_angtol.value")
    def spglib_summary(self):
        """Call spglib to find space group symmetries and Wyckoff positions."""
        s = self.structure.spget_summary(symprec=self.spglib_symprec.value,
                                         angle_tolerance=self.spglib_angtol.value)
        return pn.Row(bkw.PreText(text=s, sizing_mode='stretch_width'))

    @pn.depends("kpath_format.value", "line_density.value")
    def get_kpath(self):
        """Generate high-symmetry k-path from input structure in the ABINIT format."""
        col = pn.Column(sizing_mode='stretch_width'); ca = col.append

        s = self.structure.get_kpath_input_string(fmt=self.kpath_format.value,
                                                  line_density=self.line_density.value)
        ca(self.html_with_clipboard_btn(f"<pre> {s} </pre>"))

        if self.plot_kpath.value:
            ca("## Brillouin zone and **k**-path:")
            kpath_pane = ply(self.structure.plotly_bz(pmg_path=True, show=False), with_divider=False)
            df_kpts = dfc(self.structure.hsym_kpoints.get_highsym_datataframe(), with_divider=False)
            ca(pn.Row(kpath_pane, df_kpts))
            ca(pn.layout.Divider())

        return col

    @depends_on_btn_click('gs_input_btn')
    def on_gs_input_btn(self):
        """Generate minimalistic input file from the input structure."""
        from abipy.abio.factories import gs_input
        from abipy.data.hgh_pseudos import HGH_TABLE

        gs_inp = gs_input(
            self.structure, HGH_TABLE, kppa=self.kppra.value, ecut=8,
            spin_mode=self.label2mode[self.spin_mode.value],
            smearing=None)

        gs_inp.pop_vars(("charge", "chksymbreak"))
        gs_inp.set_vars(ecut="?? # depends on pseudos",
                        nband="?? # depends on pseudos",
                        pseudos='"pseudo1, pseudo2, ..."'
                       )

        if self.gs_type.value == "relax":
            gs_inp.set_vars(optcell=2, ionmov=2, ecutsm=0.5, dilatmx=1.05)

        gs_inp.set_mnemonics(False)

        return self.html_with_clipboard_btn(gs_inp._repr_html_())

    @depends_on_btn_click('mp_match_btn')
    def on_mp_match_btn(self):
        from abipy.core.structure import mp_match_structure
        mp = mp_match_structure(self.structure, api_key=None, endpoint=None, final=True)
        if not mp.structures:
            raise RuntimeError("No structure found in the MP database")

        return pn.Column(dfc(mp.lattice_dataframe), sizing_mode='stretch_width')

    #@depends_on_btn_click('mp_search_btn')
    #def on_mp_search_btn(self):
    #    from abipy.core.structure import mp_search
    #    chemsys_formula_id = self.stucture.formula
    #    mp = mp_search(chemsys_formula_id, api_key=None, endpoint=None, final=True)
    #    if not mp.structures:
    #        raise RuntimeError("No structure found in MP database")

    #    return pn.Column(dfc(mp.lattice_dataframe), sizing_mode='stretch_width')

    def get_panel(self, as_dict=False, **kwargs):
        """Build panel with widgets to interact with the structure either in a notebook or in a bokeh app"""
        d = {}

        d["Summary"] = pn.Row(bkw.PreText(text=self.structure.to_string(verbose=self.verbose),
                              sizing_mode="scale_both"))
        d["Spglib"] = pn.Row(
            self.pws_col(['### Spglib options', "spglib_symprec", "spglib_angtol", self.helpc("spglib_summary")]),
            self.spglib_summary
        )
        d["Kpath"] = pn.Row(
            self.pws_col(['### K-path options', "kpath_format", "line_density", "plot_kpath", self.helpc("get_kpath")]),
            self.get_kpath
        )
        d["Convert"] = pn.Row(
            self.pws_col(["### Convert structure", "output_format", self.helpc("convert")]),
            self.convert
        )
        d["Structure"] = self.get_struct_view_tab_entry()
        d["GS-input"] = pn.Row(
            self.pws_col(['### Generate GS input', "gs_type", "spin_mode", "kppra", "gs_input_btn",
                           self.helpc("on_gs_input_btn")]),
            self.on_gs_input_btn
        )
        d["MP-match"] = pn.Row(pn.Column(self.mp_match_btn), self.on_mp_match_btn)

        if as_dict: return d

        tabs = pn.Tabs(*d.items())
        return self.get_template_from_tabs(tabs, template=kwargs.get("template", None))
