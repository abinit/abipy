""""GUIs for structure."""
import param
import panel as pn
import bokeh.models.widgets as bw

from abipy.panels.core import AbipyParameterized


def _df(df, disabled=True):
    return pn.widgets.DataFrame(df, disabled=disabled)


class StructurePanel(AbipyParameterized):

    output_format = param.ObjectSelector(default="abinit",
                                         objects="abinit,cif,xsf,poscar,qe,siesta,wannier90,cssr,json".split(","),
                                         doc="Output format")

    spglib_symprec = param.Number(0.01, bounds=(0.0, None))
    spglib_angtol = param.Number(5, bounds=(0.0, None))

    kpath_format = param.ObjectSelector(default="abinit",
                                        objects=["abinit", "siesta", "wannier90"],
                                        doc="Output format")

    line_density = pn.widgets.Spinner(name="Line density", value=10, step=5, start=0, end=None)

    viewer = param.ObjectSelector(default="vesta",
                                  objects="vesta,xcrysden".split(","),
                                  doc="Viewer")
    viewer_btn = pn.widgets.Button(name="View Structure", button_type='primary')

    mp_match_btn = pn.widgets.Button(name="Connect to Materials Project", button_type='primary')

    def __init__(self, structure,  **params):
        super().__init__(**params)
        self.structure = structure

    @param.depends("output_format")
    def convert(self):
        return pn.Row(bw.PreText(text=self.structure.convert(fmt=self.output_format)), sizing_mode='stretch_width')

    @param.depends("spglib_symprec", "spglib_angtol")
    def spglib_summary(self):
        s = self.structure.spget_summary(symprec=self.spglib_symprec, angle_tolerance=self.spglib_angtol)
        return pn.Row(bw.PreText(text=s, sizing_mode='stretch_width'))

    @param.depends("kpath_format", "line_density.value")
    def get_kpath(self):
        s = self.structure.get_kpath_input_string(fmt=self.kpath_format, line_density=self.line_density.value)
        return pn.Row(bw.PreText(text=s, sizing_mode='stretch_width'))

    @param.depends("viewer_btn.clicks")
    def view(self):
        if self.viewer_btn.clicks == 0: return
        #import nglview as nv
        #view = nv.show_pymatgen(self.structure)
        #print(view)
        #print(view._display_image())
        #return view.display(gui=True)
        #return pn.interact(view)
        #return view.render_image()

        self.structure.visualize(appname=self.viewer) #.value)

    @param.depends("mp_match_btn.clicks")
    def on_mp_match_btn(self):
        if self.mp_match_btn.clicks == 0: return
        from abipy.core.structure import mp_match_structure
        mp = mp_match_structure(self.structure, api_key=None, endpoint=None, final=True)
        if not mp.structures:
            raise RuntimeError("No structure found in MP database")

        return pn.Column(_df(mp.lattice_dataframe), sizing_mode='stretch_width')

    def get_panel(self):
        """Build panel with widgets to interact with the structure either in a notebook or in a bokeh app"""
        tabs = pn.Tabs()
        w = pn.WidgetBox('# Spglib options', self.param.spglib_symprec, self.param.spglib_angtol)
        tabs.append(("Spglib", pn.Row(w, self.spglib_summary)))
        w = pn.WidgetBox('# K-path options', self.param.kpath_format, self.line_density)
        tabs.append(("Kpath", pn.Row(w, self.get_kpath)))
        tabs.append(("Convert", pn.Row(pn.Column(self.param.output_format), self.convert)))
        tabs.append(("View", pn.Row(pn.Column(self.param.viewer, self.viewer_btn), self.view)))
        tabs.append(("MP-match", pn.Row(pn.Column(self.mp_match_btn), self.on_mp_match_btn)))

        return tabs
