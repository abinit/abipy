""""GUIs for structure."""
import param
import panel as pn
import bokeh.models.widgets as bw


class StructurePanel(param.Parameterized):

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
        import nglview as nv
        view = nv.show_pymatgen(self.structure)

        #self.structure.visualize(appname="vesta")
        #print(view)
        #print(view._display_image())
        #return view.display(gui=True)
        #return pn.interact(view)
        #return view.render_image()

    def get_panel(self):
        """Build panel with widgets to interact with the structure either in a notebook or in a bokeh app"""
        tabs = pn.Tabs()
        w = pn.WidgetBox('# Spglib options', self.param.spglib_symprec, self.param.spglib_angtol)
        tabs.append(("Spglib", pn.Row(w, self.spglib_summary)))
        w = pn.WidgetBox('# K-path options', self.param.kpath_format, self.line_density)
        tabs.append(("Kpath", pn.Row(w, self.get_kpath)))
        tabs.append(("Convert", pn.Row(pn.Column(self.param.output_format), self.convert)))
        tabs.append(("View", pn.Row(pn.Column(self.param.viewer, self.viewer_btn), self.view)))
        return tabs
