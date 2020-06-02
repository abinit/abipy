""""GUIs for structure."""

import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.panels.core import AbipyParameterized, gen_id

pn.config.js_files.update({
    '$': 'https://code.jquery.com/jquery-3.4.1.slim.min.js',
    "clipboard": "https://cdn.jsdelivr.net/npm/clipboard@2/dist/clipboard.min.js",
})

pn.config.js_files["ngl"] = "https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js"


_html_with_copy_to_clipboard_ncalls = 0

def html_with_copy_to_clipboard(html, btn_cls="btn-primary btn-sm"):
    global _html_with_copy_to_clipboard_ncalls
    _html_with_copy_to_clipboard_ncalls += 1
    myid = gen_id()
    html = f"""
<button class="btn {btn_cls}" data-clipboard-target="#{myid}"> Copy to clipboard </button>
<div id="{myid}"> {html} </div>
"""
    if _html_with_copy_to_clipboard_ncalls == 1:
        html += "<script>$(document).ready(function () {new ClipboardJS('.btn')})</script>"

    return pn.pane.HTML(html)


class StructurePanel(AbipyParameterized):
    """
    Panel with widgets to interact with an AbiPy Structure
    """
    verbose = param.Integer(0, bounds=(0, None), doc="Verbosity Level")

    # Convert widgets.
    output_format = pnw.Select(name="format", value="abinit",
                               options="abinit,cif,xsf,poscar,qe,siesta,wannier90,cssr,json".split(","))

    # Spglib widgets
    spglib_symprec = pnw.Spinner(name="symprec", value=0.01, start=0.0, end=None, step=0.01)
    spglib_angtol = pnw.Spinner(name="angtol", value=5, start=0.0, end=None, step=1)

    # K-path widgets
    kpath_format = pnw.Select(name="format", value="abinit", options=["abinit", "siesta", "wannier90"])
    line_density = pnw.Spinner(name="line density", value=10, step=5, start=0, end=None)

    # Viewer widgets.
    viewer_btn = pnw.Button(name="View structure", button_type='primary')
    viewer = pnw.Select(name="Viewer", value="vesta",
                        options=["vesta", "xcrysden", "vtk", "matplotlib", "mayavi"])

    # Mp-match
    mp_match_btn = pnw.Button(name="Connect to Materials Project", button_type='primary')

    # Mp-search
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
        self.structure = structure

    @param.depends("output_format.value")
    def convert(self):
        return pn.Row(bkw.PreText(text=self.structure.convert(fmt=self.output_format.value)),
                      sizing_mode='stretch_width')

    @param.depends("spglib_symprec.value", "spglib_angtol.value")
    def spglib_summary(self):
        s = self.structure.spget_summary(symprec=self.spglib_symprec.value,
                                         angle_tolerance=self.spglib_angtol.value)
        return pn.Row(bkw.PreText(text=s, sizing_mode='stretch_width'))

    @param.depends("kpath_format.value", "line_density.value")
    def get_kpath(self):
        s = self.structure.get_kpath_input_string(fmt=self.kpath_format.value,
                                                  line_density=self.line_density.value)
        return pn.Row(bkw.PreText(text=s, sizing_mode='stretch_width'))

    @param.depends("viewer_btn.clicks")
    def view(self):
        if self.viewer_btn.clicks == 0: return

        import panel as pn

        #view = self.structure.get_jsmol_view()
        #from ipywidgets_bokeh import IPyWidget
        #view = IPyWidget(widget=view) #, width=800, height=300)
        #from IPython.display import display
        #display(view)
        #return pn.Row(display(view))
        #return pn.panel(view)

        js_files = {'ngl': 'https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js'}
        pn.extension(comms='ipywidgets', js_files=js_files)
        view = self.structure.get_ngl_view()
        return pn.panel(view)

        #pn.config.js_files["ngl"]="https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js"
        #pn.extension()

        html = """<div id="viewport" style="width:100%; height:100%;"></div>
        <script>
        stage = new NGL.Stage("viewport");
        stage.loadFile("rcsb://1NKT.mmtf", {defaultRepresentation: true});
        </script>"""

#        html = """
#         <script>
#    document.addeventlistener("domcontentloaded", function () {
#      var stage = new ngl.stage("viewport");
#      stage.loadfile("rcsb://1crn", {defaultrepresentation: true});
#    });
#  </script>"""

#        html = """
#<script>
#document.addeventlistener("domcontentloaded", function () {
#    // create a `stage` object
#    var stage = new NGL.Stage("viewport");
#    // load a PDB structure and consume the returned `Promise`
#    stage.loadFile("rcsb://1CRN").then(function (component) {
#    // add a "cartoon" representation to the structure component
#    component.addRepresentation("cartoon");
#    // provide a "good" view of the structure
#    component.autoView();
#  });
#});
#</script>"""

        ngl_pane = pn.pane.HTML(html, height=500, width=500)
        return pn.Row(ngl_pane)

        view = self.structure.get_ngl_view()
        #return self.structure.crystaltoolkitview()
        #import nglview as nv
        #view = nv.demo(gui=False)


        #from bokeh.models import ColumnDataSource
        #from bokeh.io import show, curdoc
        #from bokeh.models.widgets import Button, TextInput
        #from bokeh.layouts import layout, widgetbox
        #from jsmol_bokeh_extension import JSMol
        #script_source = ColumnDataSource()

        #info = dict(
        #    height="100%",
        #    width="100%",
        #    serverURL="https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
        #    use="HTML5",
        #    j2sPath="https://chemapps.stolaf.edu/jmol/jsmol/j2s",
        #    script=
        #    "background black;load https://chemapps.stolaf.edu/jmol/jsmol-2013-09-18/data/caffeine.mol",
        #)

        #applet = JSMol(
        #    width=600,
        #    height=600,
        #    script_source=script_source,
        #    info=info,
        #)

        #button = Button(label='Execute')
        #inp_script = TextInput(value='background white;')

        #def run_script():
        #    script_source.data['script'] = [inp_script.value]

        #button.on_click(run_script)
        #ly = layout([applet, widgetbox(button, inp_script)])
        #show(ly)
        #curdoc().add_root(ly)

        return pn.Row(applet)
        return self.structure.visualize(appname=self.viewer.value)

    @param.depends("gs_input_btn.clicks")
    def on_gs_input_btn(self):
        if self.gs_input_btn.clicks == 0: return
        from abipy.abio.factories import gs_input
        from abipy.data.hgh_pseudos import HGH_TABLE

        gs_inp = gs_input(
            self.structure, HGH_TABLE, kppa=self.kppra.value, ecut=8,
            spin_mode=self.label2mode[self.spin_mode.value],
            smearing=None)
        gs_inp.pop_vars(("charge", "chksymbreak"))
        gs_inp.set_vars(ecut="?? # depends on pseudos", nband="?? # depends on pseudos")

        if self.gs_type.value == "relax":
            gs_inp.set_vars(optcell=2, ionmov=2, ecutsm=0.5, dilatmx=1.05)

        gs_inp.set_mnemonics(False)
        return html_with_copy_to_clipboard(gs_inp._repr_html_())

    @param.depends("mp_match_btn.clicks")
    def on_mp_match_btn(self):
        if self.mp_match_btn.clicks == 0: return
        from abipy.core.structure import mp_match_structure
        mp = mp_match_structure(self.structure, api_key=None, endpoint=None, final=True)
        if not mp.structures:
            raise RuntimeError("No structure found in MP database")

        return pn.Column(self._df(mp.lattice_dataframe), sizing_mode='stretch_width')

    #@param.depends("mp_search_btn.clicks")
    #def on_mp_search_btn(self):
    #    if self.mp_search_btn.clicks == 0: return
    #    from abipy.core.structure import mp_search
    #    chemsys_formula_id = self.stucture.formula
    #    mp = mp_search(chemsys_formula_id, api_key=None, endpoint=None, final=True)
    #    if not mp.structures:
    #        raise RuntimeError("No structure found in MP database")

    #    return pn.Column(self._df(mp.lattice_dataframe), sizing_mode='stretch_width')

    def get_panel(self):
        """Build panel with widgets to interact with the structure either in a notebook or in a bokeh app"""
        tabs = pn.Tabs(); app = tabs.append
        row = pn.Row(bkw.PreText(text=self.structure.to_string(verbose=self.verbose), sizing_mode="scale_both"))
        app(("Summary", row))
        ws = pn.Column('# Spglib options', self.spglib_symprec, self.spglib_angtol)
        app(("Spglib", pn.Row(ws, self.spglib_summary)))
        ws = pn.Column('# K-path options', self.kpath_format, self.line_density)
        app(("Kpath", pn.Row(ws, self.get_kpath)))
        app(("Convert", pn.Row(pn.Column("# Convert structure", self.output_format), self.convert)))
        app(("View", pn.Row(pn.Column("# Visualize structure", self.viewer, self.viewer_btn), self.view)))
        ws = pn.Column('# Generate GS input', self.gs_type, self.spin_mode, self.kppra, self.gs_input_btn)
        app(("GS-input", pn.Row(ws, self.on_gs_input_btn)))
        app(("MP-match", pn.Row(pn.Column(self.mp_match_btn), self.on_mp_match_btn)))

        return tabs
