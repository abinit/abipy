""""GUIs for structure."""

import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.panels.core import HasStructureParams, ButtonContext, dfc #gen_id,

#pn.config.js_files.update({
#    #'$': 'https://code.jquery.com/jquery-3.4.1.slim.min.js',
#    "clipboard": "https://cdn.jsdelivr.net/npm/clipboard@2/dist/clipboard.min.js",
#})

#pn.config.js_files["ngl"] = "https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js"


class StructurePanel(HasStructureParams):
    """
    Panel with widgets to interact with an AbiPy Structure
    """
    # Convert widgets.
    output_format = pnw.Select(name="format", value="abinit",
                               options="abinit,cif,xsf,poscar,qe,siesta,wannier90,cssr,json".split(","))

    # Spglib widgets
    spglib_symprec = pnw.Spinner(name="symprec", value=0.01, start=0.0, end=None, step=0.01)
    spglib_angtol = pnw.Spinner(name="angtol", value=5, start=0.0, end=None, step=1)

    # K-path widgets
    kpath_format = pnw.Select(name="format", value="abinit", options=["abinit", "siesta", "wannier90"])
    line_density = pnw.Spinner(name="line density", value=10, step=5, start=0, end=None)

    # MP-match
    mp_match_btn = pnw.Button(name="Connect to Materials Project", button_type='primary')

    # MP-search
    #mp_search_btn = pnw.Button(name="Connect to Materials Project", button_type='primary')
    #mp_api_key

    # GS input generator widgets.
    gs_input_btn = pnw.Button(name="Generate input", button_type='primary')
    gs_type = pnw.Select(name="GS type", value="scf", options=["scf", "relax"])
    kppra = pnw.Spinner(name="kppra", value=1000, step=1000, start=0, end=None)

    label2mode = {
        "unpolarized": 'unpolarized',
        "polarized": 'polarized',
        "anti-ferromagnetic": "afm",
        "non-collinear with magnetism": "spinor",
        "non-collinear, no magnetism": "spinor_nomag",
    }

    spin_mode = pnw.Select(name="SpinMode", value="unpolarized", options=list(label2mode.keys()))

    def __init__(self, structure,  **params):
        super().__init__(**params)
        self._structure = structure

    @property
    def structure(self):
        return self._structure

    @param.depends("output_format.value")
    def convert(self):
        """Convert the input structure to one of the format selected by the user."""
        s = self.structure.convert(fmt=self.output_format.value)
        return self.html_with_clipboard_btn(f"<pre> {s} </pre>")

    @param.depends("spglib_symprec.value", "spglib_angtol.value")
    def spglib_summary(self):
        """Call spglib to find space group symmetries and Wyckoff positions."""
        s = self.structure.spget_summary(symprec=self.spglib_symprec.value,
                                         angle_tolerance=self.spglib_angtol.value)
        return pn.Row(bkw.PreText(text=s, sizing_mode='stretch_width'))

    @param.depends("kpath_format.value", "line_density.value")
    def get_kpath(self):
        """Generate high-symmetry k-path from input structure in ABINIT format.."""
        s = self.structure.get_kpath_input_string(fmt=self.kpath_format.value,
                                                  line_density=self.line_density.value)
        return self.html_with_clipboard_btn(f"<pre> {s} </pre>")

    @param.depends("gs_input_btn.clicks")
    def on_gs_input_btn(self):
        """Generate minimalistic input file from the input structure."""
        if self.gs_input_btn.clicks == 0: return

        with ButtonContext(self.gs_input_btn):
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

    @param.depends("mp_match_btn.clicks")
    def on_mp_match_btn(self):
        """Match input structure with the structures available on the materials project."""
        if self.mp_match_btn.clicks == 0: return

        with ButtonContext(self.mp_match_btn):
            from abipy.core.structure import mp_match_structure
            mp = mp_match_structure(self.structure, api_key=None, endpoint=None, final=True)
            if not mp.structures:
                raise RuntimeError("No structure found in the MP database")

            return pn.Column(dfc(mp.lattice_dataframe), sizing_mode='stretch_width')

    #@param.depends("mp_search_btn.clicks")
    #def on_mp_search_btn(self):
    #    if self.mp_search_btn.clicks == 0: return
    #    from abipy.core.structure import mp_search
    #    chemsys_formula_id = self.stucture.formula
    #    mp = mp_search(chemsys_formula_id, api_key=None, endpoint=None, final=True)
    #    if not mp.structures:
    #        raise RuntimeError("No structure found in MP database")

    #    return pn.Column(dfc(mp.lattice_dataframe), sizing_mode='stretch_width')

    def get_panel(self):
        """Build panel with widgets to interact with the structure either in a notebook or in a bokeh app"""

        tabs = pn.Tabs(); app = tabs.append

        app(("Summary",
            pn.Row(bkw.PreText(text=self.structure.to_string(verbose=self.verbose), sizing_mode="scale_both"))
        ))
        app(("Spglib", pn.Row(
            pn.Column('# Spglib options', *self.pws("spglib_symprec", "spglib_angtol", self.helpc("spglib_summary"))),
            self.spglib_summary)
        ))
        app(("Kpath", pn.Row(
            pn.Column('# K-path options', *self.pws("kpath_format", "line_density", self.helpc("get_kpath"))),
            self.get_kpath)
        ))
        app(("Convert", pn.Row(
            pn.Column("# Convert structure", *self.pws("output_format", self.helpc("convert"))),
            self.convert)
        ))
        app(self.get_struct_view_tab_entry())
        app(("GS-input", pn.Row(
            pn.Column('# Generate GS input', *self.pws("gs_type", "spin_mode", "kppra", "gs_input_btn",
                      self.helpc("on_gs_input_btn"))),
            self.on_gs_input_btn)
        ))
        app(("MP-match", pn.Row(
            pn.Column(self.mp_match_btn),
            self.on_mp_match_btn)
        ))

        return tabs
