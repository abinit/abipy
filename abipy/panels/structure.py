""""GUIs for structure."""

import os
import io
import param
import panel as pn
import numpy as np
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from abipy.core.structure import Structure
from abipy.panels.core import AbipyParameterized, HasStructureParams, dfc, mpl, ply, depends_on_btn_click, Loading


def _make_targz_bytes(inp_or_multi, remove_dir=True):
    """
    Return bytesIO with the targz file containing input and pseudos.
    """
    targz_path = inp_or_multi.make_targz()
    output = io.BytesIO()
    with open(targz_path, "rb") as fh:
        output.write(fh.read())
    output.seek(0)

    # TODO Remove directory. Implement with in make_targz.
    #if remove_dir:

    return output


class StructurePanel(HasStructureParams):
    """
    Panel with widgets to interact with an AbiPy Structure
    """

    def __init__(self, structure,  with_inputs=True, **params):
        """
        Args:
            structure: |Structure| object.
            with_inputs: True if tabs for generating input files should be shown.
        """
        self._structure = structure
        self.with_inputs = with_inputs

        # Convert widgets.
        self.output_format = pnw.Select(name="format", value="abinit",
                                        options="abinit,cif,xsf,poscar,qe,siesta,wannier90,cssr,json".split(","))

        # Spglib widgets
        self.spglib_symprec = pnw.Spinner(name="symprec", value=0.01, start=0.0, end=None, step=0.01)
        self.spglib_angtol = pnw.Spinner(name="angtol", value=5, start=0.0, end=None, step=1)

        # Abisanitize widgets
        self.abisanitize_btn = pnw.Button(name="Run abisanitize", button_type='primary')
        self.select_primitive = pnw.Select(name='Select primitive',
                                           options=['primitive', 'primitive_standard', "no_primitive"])

        # K-path widgets
        self.kpath_format = pnw.Select(name="format", value="abinit", options=["abinit", "siesta", "wannier90"])
        self.line_density = pnw.Spinner(name="line density", value=10, step=5, start=0, end=None)
        self.plot_kpath = pnw.Checkbox(name='Plot k-path', value=False)

        # MP-match
        self.mp_match_btn = pnw.Button(name="Connect to Materials Project", button_type='primary')

        # MP-search
        #mp_search_btn = pnw.Button(name="Connect to Materials Project", button_type='primary')
        #mp_api_key

        # Widgets to control the generation of input files.
        self.kppra = pnw.Spinner(name="kppra", value=1000, step=500, start=0, end=None)
        self.smearing_type = pnw.Select(name="Smearing type", value=None,
                                        options=[None, "gaussian", "fermi_dirac"])
        self.tsmear = pnw.Spinner(name="tsmear (Ha)", value=0.01, step=0.002, start=0.0, end=None)

        self.xc_type = pnw.Select(name="XC type", value="LDA", options=["LDA", "PBEsol", "PBE"])
        self.pseudos_type = pnw.Select(name="Pseudos type", value="NC", options=["NC", "PAW"])

        # widgets for GS input generator
        self.gs_type = pnw.Select(name="GS type", value="scf", options=["scf", "relax"])
        self.gs_input_btn = pnw.Button(name="Generate input", button_type='primary')

        # widgets for e-bands input generator.
        self.ebands_input_btn = pnw.Button(name="Generate input", button_type='primary')
        self.edos_kppra = pnw.Spinner(name="edos_kppra", value=0, step=500, start=0, end=None)

        # widgets for DFPT phonons input generator.
        self.ph_input_btn = pnw.Button(name="Generate input", button_type='primary')
        self.with_becs = pnw.Checkbox(name='Born effective charges (polar semiconductor)', value=False)
        #self.dfpt_ngqpt = pnw.LiteralInput(name='ngqpt (python list)', value=[None, None, None], type=list)

        self.label2mode = {
            "unpolarized": 'unpolarized',
            "polarized": 'polarized',
            "non-collinear SOC with magnetism": "spinor",
            "non-collinear SOC, no magnetism": "spinor_nomag",
            "collinear anti-ferromagnetic": "afm",
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

    @depends_on_btn_click('abisanitize_btn')
    def on_abisanitize_btn(self):
        """
        Returns a new structure in which:

            - Structure is refined.
            - Reduced to primitive settings.
            - Lattice vectors are exchanged if the triple product is negative

        **symprec**:
            Symmetry precision used to refine the structure.

        **angle_tolerance**:
            Tolerance on angles.

        **primitive**:
            Returns most primitive structure found.

        **primitive_standard**:
            Whether to convert to a primitive cell using the standards defined in Setyawan, W., & Curtarolo, S. (2010).
            High-throughput electronic band structure calculations: Challenges and tools.
            Computational Materials Science, 49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
        """
        primitive, primitive_standard = False, False
        if self.select_primitive.value in ('primitive', 'primitive_standard'):
            primitive = True
            primitive_standard = self.select_primitive.value == 'primitive_standard'
        else:
            assert self.select_primitive.value == "no_primitive"

        #print("primitive", primitive, "primitive_standard", primitive_standard)

        s = self.structure.abi_sanitize(symprec=self.spglib_symprec.value,
                                        angle_tolerance=self.spglib_angtol.value,
                                        primitive=primitive, primitive_standard=primitive_standard)

        return pn.Row(
                pn.Column("## Input Structure:", bkw.PreText(text=str(self.structure), sizing_mode='stretch_width')),
                pn.Column("## Sanitized:", bkw.PreText(text=str(s), sizing_mode='stretch_width')),
                sizing_mode='stretch_width')

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

    def _get_pseudos(self):
        #self.xc_type.value
        #self.pseudos_type.value

        from abipy.data.hgh_pseudos import HGH_TABLE
        return HGH_TABLE

        #if self.pseudo_type.value == "NC":
        #elif self.pseudo_type.value == "PAW":
        #    raise ValueError("PAW pseudos are not available!")

    def _get_smearing(self):
        smearing = None
        if self.smearing_type.value is not None:
            smearing = "%s:%s Ha" % (self.smearing_type.value, self.tsmear.value)
        return smearing

    def _finalize(self, inp_or_multi):

        # TODO ckksymbreak should be set to 0 else Abinit may stop.
        # DEBUG MODE
        #if inp_or_multi.is_input:
        #    this = inp_or_multi.deepcopy()
        #    this.pop_vars("pseudos")
        #    v = this.abivalidate()
        #    v_list = [v]
        #else:
        #    this = inp_or_multi.deepcopy()
        #    this.pop_vars("pseudos")
        #    v_list = this.abivalidate()

        #count = 0
        #for v in v_list:
        #    if v.retcode != 0:
        #        count += 1
        #        print(v.stderr_file.read())
        #        print(v.log_file.read())
        #if count != 0: raise RuntimeError("")


        # Here we replace the absolute paths with their basenames
        # so that we can create a targz file with the pseudos
        # We do it at this level because absolute paths are needed for executing abinit.
        inp_or_multi.set_vars(pseudos='"%s"' % ", ".join(p.basename for p in inp_or_multi.pseudos))

        html = self.html_with_clipboard_btn(inp_or_multi._repr_html_())
        #return html

        def download_input():
            return _make_targz_bytes(inp_or_multi)

        file_download = pnw.FileDownload(filename="input.tar.gz", callback=download_input)

        msg = """
This input file has been automatically generated by AbiPy.
Several input parameters have **default values** that might not be suitable for you particular calculation.
Please check the input variable and modify it according to your needs.

Also note that running big calculations with lot of datasets is not the most efficient approach.
Examples of AbiPy scripst to automate calculations without datasets are available
[on this page](https://abinit.github.io/abipy/flow_gallery/index.html)
"""
        alert = pn.pane.Alert(msg, alert_type="danger")

        return pn.Column(html, file_download, pn.layout.Divider(), alert,
                         sizing_mode="stretch_width")

    def get_gs_input(self):
        """
        Return an AbinitInput for GS calculation from the parameters selected via the widgets.
        """
        from abipy.abio.factories import gs_input
        gs_inp = gs_input(structure=self.structure, pseudos=self._get_pseudos(),
                          kppa=self.kppra.value,
                          ecut=8,
                          spin_mode=self.label2mode[self.spin_mode.value],
                          smearing=self._get_smearing(),
                          charge=0.0,
                         )

        gs_inp.set_mnemonics(False)

        # TODO: Should reorder pseudos?
        gs_inp.pop_vars(("charge", "chksymbreak"))
        #gs_inp.set_vars(#ecut="??  # depends on pseudos",
        #                #nband="?? # depends on pseudos",
        #                pseudos='"%s"' % ", ".join(p.basename for p in gs_inp.pseudos),
        #                )

        return gs_inp

    @depends_on_btn_click('gs_input_btn')
    def on_gs_input_btn(self):
        """
        Generate minimalistic input file from the structure.
        """
        gs_inp = self.get_gs_input()

        if self.gs_type.value == "relax":
            gs_inp.set_vars(optcell=2, ionmov=2, ecutsm=0.5, dilatmx=1.05)

        return self._finalize(gs_inp)

    @depends_on_btn_click('ebands_input_btn')
    def on_ebands_input_btn(self):
        """
        Generate minimalistic input file from band structure calculation.
        """
        from abipy.abio.factories import ebands_input

        dos_kppa = self.edos_kppra
        if dos_kppa == 0.0: dos_kppa = None

        multi = ebands_input(structure=self.structure, pseudos=self._get_pseudos(),
                             kppa=self.kppra.value,
                             nscf_nband=None,
                             ndivsm=10,
                             ecut=8,
                             pawecutdg=None,
                             scf_nband=None,
                             spin_mode=self.label2mode[self.spin_mode.value],
                             smearing=self._get_smearing(),
                             charge=0.0,
                             dos_kppa=None,
        )

        # Add getwfk variables.
        for inp in multi[1:]:
            inp["getwfk"] = 1

        multi.pop_vars(("charge", "chksymbreak"))
        #multi.set_vars(#ecut="??  # depends on pseudos",
        #               #nband="?? # depends on pseudos",
        #               pseudos='"%s"' % ", ".join(p.basename for p in multi.pseudos),
        #               )

        return self._finalize(multi)

    @depends_on_btn_click('ph_input_btn')
    def on_ph_input_btn(self):
        """
        Generate minimalistic input file from phonon calculations with DFPT.
        """
        #raise ValueError("hello foobar")
        from abipy.abio.factories import phonons_from_gsinput

        gs_inp = self.get_gs_input()

        ph_ngqpt = [2, 2, 2]
        #ph_ngqpt = gs_inp["ngkpt"]
        qpoints = gs_inp.abiget_ibz(ngkpt=ph_ngqpt, shiftk=[0, 0, 0], kptopt=1).points

        ndtset = 1 + len(qpoints)
        qstart = 1
        with_becs = self.with_becs.value
        if with_becs:
            # Account for the DDK (2nd dataset)
            ndtset += 1
            qstart = 2

        multi = gs_inp.replicate(ndtset)

        if with_becs:
            ddk_input = multi[qstart - 1]
            ddk_input.set_vars(
                    rfelfd=2,          # only the derivative of ground-state wavefunctions with respect to k
                    rfdir=[1, 1, 1],
                    nqpt=1,
                    qpt=(0, 0, 0),
                    kptopt=2,         # 2 to take into account time-reversal symmetry.
                    iscf=-3,          # The d/dk perturbation must be treated in a non-self-consistent way
                    getwfk=1,
                    comment="Response Function calculation: d/dk",
            )
            ddk_input.pop_tolerances()
            ddk_input.set_vars(tolwfr=1.0e-22)

        for qpt, ph_input in zip(qpoints, multi[qstart:]):
            is_gamma = np.sum(qpt ** 2) < 1e-12
            #if with_becs and is_gamma: continue
            ph_input.set_vars(
                    kptopt=2 if is_gamma else 3, # use time-reversal if Gamma
                    rfphon=1,                    # Will consider phonon-type perturbation
                    nqpt=1,                      # One wavevector is to be considered
                    qpt=qpt,                     # q-wavevector.
                    rfatpol=[1, len(gs_inp.structure)],
                    rfdir=[1, 1, 1],
                    getwfk=1,
                    rfelfd=3 if (with_becs and is_gamma) else None,
                    getddk=2 if (with_becs and is_gamma) else None,
                    prtwf=-1,
                    comment="Input file for PH calculation with DFPT.",
            )
            ph_input.pop_tolerances()
            ph_input.set_vars(tolvrs=1.0e-10)


        #multi = phonons_from_gsinput(gs_inp, ph_ngqpt=None, qpoints=None, with_ddk=True,
        #                             with_dde=True, with_bec=False,
        #                             ph_tol=None, ddk_tol=None, dde_tol=None, wfq_tol=None,
        #                             qpoints_to_skip=None, manager=None)

        # TODO: Should add gs_inp to multi
        # Add getwfk variables.
        #for inp in multi[1:]:
        #    inp["getwfk"] = 1

        #multi.pop_vars(("charge", "chksymbreak"))
        #multi.set_vars(#ecut="??  # depends on pseudos",
        #               #nband="?? # depends on pseudos",
        #               pseudos='"%s"' % ", ".join(p.basename for p in multi.pseudos),
        #               )

        return self._finalize(multi)

    @depends_on_btn_click('mp_match_btn')
    def on_mp_match_btn(self):
        mp = self.structure.mp_match()
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
        d["AbiSanitize"] = pn.Row(
            self.pws_col(['### Spglib options', "spglib_symprec", "spglib_angtol", "select_primitive", "abisanitize_btn",
                          pn.layout.Divider(), self.helpc("on_abisanitize_btn")]),
            self.on_abisanitize_btn
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
            self.pws_col(['### Generate GS input', "gs_type",
                          "spin_mode", "kppra", "smearing_type", "tsmear",
                          "xc_type", "pseudos_type",
                          "gs_input_btn",
                           self.helpc("on_gs_input_btn")]),
            self.on_gs_input_btn
        )

        # Add tabs to generate inputs from structure.
        if self.with_inputs:
            d["Ebands-input"] = pn.Row(
                self.pws_col(['### Generate Ebands input',
                              "spin_mode", "kppra", "edos_kppra", "smearing_type", "tsmear",
                              "xc_type", "pseudos_type",
                              "ebands_input_btn",
                               self.helpc("on_ebands_input_btn")]),
                self.on_ebands_input_btn
            )
            d["PH-input"] = pn.Row(
                self.pws_col(['### Generate phonon input',
                              "spin_mode", "kppra", "smearing_type", "tsmear",
                              "with_becs",
                              "xc_type", "pseudos_type",
                              "ph_input_btn",
                               self.helpc("on_ph_input_btn")]),
                self.on_ph_input_btn
            )
            d["MP-match"] = pn.Column(pn.Row(self.mp_match_btn, align="center"),
                                      pn.layout.Divider(),
                                      self.on_mp_match_btn,
                                      sizing_mode="stretch_width")

        if as_dict: return d

        tabs = pn.Tabs(*d.items())
        return self.get_template_from_tabs(tabs, template=kwargs.get("template", None))


class InputFileGenerator(AbipyParameterized):

    #symprec
    #angle_tolerance

    def __init__(self, **params):

        super().__init__(**params)

        # Spglib widgets
        #self.spglib_symprec = pnw.Spinner(name="symprec", value=0.01, start=0.0, end=None, step=0.01)
        #self.spglib_angtol = pnw.Spinner(name="angtol", value=5, start=0.0, end=None, step=1)

        help_md = pn.pane.Markdown("""
## Main area

This web app exposes some of the post-processing capabilities of AbiPy.
Use the **Choose File** to upload one of the files supported by this app.
Keep in mind that **the file extension matters!**
""")

        self.main_area = pn.Column(help_md,
                                   self.get_alert_data_transfer(),
                                   sizing_mode="stretch_width")

        self.file_input = pnw.FileInput(height=60, css_classes=["pnx-file-upload-area"])
        self.file_input.param.watch(self.on_file_input, "value")

        self.mpid_input = pnw.TextInput(name='mp-id', placeholder='Enter e.g. mp-149 for Silicon and press ⏎')
        self.mpid_input.param.watch(self.on_mpid_input, "value")
        self.mpid_input_errstr = pn.Str("")

    def on_file_input(self, event):
        print("filename", self.file_input.filename)
        if self.file_input.value is None: return None
        #print("value", self.file_input.value)
        import tempfile
        workdir = tempfile.mkdtemp()

        fd, tmp_path = tempfile.mkstemp(suffix=self.file_input.filename)
        with open(tmp_path, "wb") as fh:
            fh.write(self.file_input.value)

        self.input_structure = Structure.from_file(tmp_path)
        os.remove(tmp_path)
        self.update_main_area()

    def on_mpid_input(self, event):
        with Loading(self.mpid_input):
            self.input_structure = Structure.from_mpid(self.mpid_input.value)

        self.update_main_area()

    def update_main_area(self):
        #self.structure = self.input_structure.abi_sanitize(symprec=1e-3, angle_tolerance=5,
        #                                                   primitive=True, primitive_standard=False)

        d = self.input_structure.get_panel(as_dict=True)
        d = {k: d[k] for k in ("GS-input", "Ebands-input", "PH-input", "Summary")}
        tabs = pn.Tabs(*d.items())
        self.main_area.objects = [tabs]
        #self.main_area.objects = [self.input_structure.get_panel()]

    def get_panel(self):

        col = pn.Column(
            "## Upload **any file** with a structure (*.nc*, *.abi*, *.cif*, *.xsf*, *POSCAR*):",
            self.get_fileinput_section(self.file_input),
            "## or get the structure from the [Materials Project](https://materialsproject.org/) database:",
            pn.Row(self.mpid_input, self.mpid_input_errstr),
            sizing_mode="stretch_width")

        main = pn.Column(col, self.main_area, sizing_mode="stretch_width")
        cls, kwds = self.get_abinit_template_cls_kwds()

        return cls(main=main, title="Input File Generator", **kwds)
