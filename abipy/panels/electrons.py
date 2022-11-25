""""
AbiPy panels for electronic properties.
"""

import param
import panel as pn
import panel.widgets as pnw

from abipy.panels.core import AbipyParameterized, ActiveBar, Loading, ply, mpl, depends_on_btn_click


class CompareEbandsWithMP(AbipyParameterized):

    with_gaps = param.Boolean(True)

    ylims_ev = param.Range(default=(-10, +10), doc="Energy window around the Fermi energy.")

    info_str = """
This app alllows users to upload two files with KS energies.
"""

    def __init__(self, **params):

        super().__init__(**params)

        help_md = pn.pane.Markdown(f"""
## Description

{self.info_str}
""")

        self.main_area = pn.Column(help_md,
                                   self.get_alert_data_transfer(),
                                   sizing_mode="stretch_width")

        self.replot_btn = pnw.Button(name="Replot", button_type='primary')

        self.file_input = pnw.FileInput(height=60, css_classes=["pnx-file-upload-area"])
        self.file_input.param.watch(self.on_file_input, "value")
        self.mp_progress = pn.indicators.Progress(name='Fetching data from the MP website', bar_color="warning",
                                                  active=False, width=200, height=10, align="center")

    def on_file_input(self, event):
        self.abinit_ebands = self.get_ebands_from_file_input(self.file_input)

        # Match Abinit structure with MP
        mp = self.abinit_ebands.structure.mp_match()
        if not mp.structures:
            raise RuntimeError("No structure found in the MP database")

        # Get structures from MP as AbiPy ElectronBands.
        from abipy.electrons.ebands import ElectronBands
        self.mp_ebands_list = []
        with ActiveBar(self.mp_progress):
            for mp_id in mp.ids:
                if mp_id == "this": continue
                eb = ElectronBands.from_mpid(mp_id)
                self.mp_ebands_list.append(eb)

        self.update_main()

    def update_main(self):
        with Loading(self.main_area):
            col = self.pws_col(["## Plot options", "with_gaps", "ylims_ev", "replot_btn"])
            ca = col.append

            ca("## Abinit Electronic band structure:")
            ylims = self.ylims_ev
            ca(ply(self.abinit_ebands.plotly(e0="fermie", ylims=ylims, with_gaps=self.with_gaps, show=False)))

            for mp_ebands in self.mp_ebands_list:
                ca("## MP Electronic band structure:")
                ca(ply(mp_ebands.plotly(e0="fermie", ylims=ylims, with_gaps=self.with_gaps, show=False)))

            #self.main_area.objects = [col]
            self.main_area.objects = col.objects

    @depends_on_btn_click('replot_btn')
    def on_replot_btn(self):
        self.update_main()

    def get_panel(self):
        col = pn.Column(
            "## Upload a *nc* file with energies along a **k**-path (possibly a *GSR.nc* file):",
            self.get_fileinput_section(self.file_input),
            pn.Row("## Fetching data from the MP website: ", self.mp_progress,
                   sizing_mode="stretch_width"),
            sizing_mode="stretch_width")

        main = pn.Column(col, self.main_area, sizing_mode="stretch_width")

        cls, kwds = self.get_abinit_template_cls_kwds()

        return cls(main=main, title="Compare with MP Ebands", **kwds)


class SkwPanelWithFileInput(AbipyParameterized):

    lpratio = param.Integer(default=5, bounds=(1, None),
                            doc="Ratio between number of k-points and number of star-functions")

    info_str = """
This app allows users to upload two files with KS energies.
The first file contains the energies in the IBZ used for the SKW interpolation (NB: this file is required).
The second (optional) file contains the energies along a k-path.
The interpolated energies are then compared with the ab-initio ones on the k-path.
The user can change the SKW intepolation parameters to gauge the quality of the SKW fit.
"""

    def __init__(self, **params):

        super().__init__(**params)

        help_md = pn.pane.Markdown(f"""
## Description

{self.info_str}
""")

        self.main_area = pn.Column(help_md,
                                   self.get_alert_data_transfer(),
                                   sizing_mode="stretch_width")

        self.ibz_file_input = pnw.FileInput(height=60, css_classes=["pnx-file-upload-area"])
        self.ibz_file_input.param.watch(self.on_ibz_file_input, "value")
        self.ebands_ibz = None

        self.kpath_file_input = pnw.FileInput(height=60, css_classes=["pnx-file-upload-area"])
        self.kpath_file_input.param.watch(self.on_kpath_file_input, "value")
        self.ebands_kpath = None

    def on_ibz_file_input(self, event):
        self.ebands_ibz = self.get_ebands_from_file_input(self.ibz_file_input)
        self.update_main_area()

    def on_kpath_file_input(self, event):
        self.ebands_kpath = self.get_ebands_from_file_input(self.kpath_file_input)
        self.update_main_area()

    def update_main_area(self):

        with Loading(self.main_area):

            if self.ebands_kpath is None or self.ebands_ibz is None: return

            # SKW interpolation
            r = self.ebands_ibz.interpolate(lpratio=self.lpratio, filter_params=None)

            # Build plotter.
            plotter = self.ebands_kpath.get_plotter_with("Ab-initio", "SKW interp", r.ebands_kpath)

            mpl_pane = mpl(plotter.combiplot(**self.mpl_kwargs))

            col = pn.Column(mpl_pane, sizing_mode="stretch_width")

            self.main_area.objects = [col]

    def get_panel(self):
        col = pn.Column(
            "## Upload (or drag & drop) any *nc* file with energies in the IBZ (_possibly a *GSR.nc* file_):",
            self.get_fileinput_section(self.ibz_file_input),
            "## Upload (or drag & drop) any *nc* file with energies along a **k**-path (_possibly a *GSR.nc* file_):",
            self.get_fileinput_section(self.kpath_file_input),
            sizing_mode="stretch_width")

        main = pn.Column(col, self.main_area, sizing_mode="stretch_width")
        cls, kwds = self.get_abinit_template_cls_kwds()

        return cls(main=main, title="SKW Analyzer", **kwds)
