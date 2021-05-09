""""Basic tools and mixin classes for AbiPy panels."""

import io
#import pathlib
import tempfile
import numpy as np
import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw
import pandas as pd

from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core import abinit_units as abu
from abipy.tools.plotting import push_to_chart_studio


#pathlib.Path(__file__) / "assets"


def abipanel():
    """
    Activate panel extensions used by AbiPy. Return panel module.
    """

    try:
        import panel as pn
    except ImportError as exc:
        cprint("Use `conda install panel` or `pip install panel` to install the python package.", "red")
        raise exc

    extensions = [
        "plotly",
        #"mathjax",
        #"katex",
    ]

    #print("loading extensions:", extensions)
    pn.extension(*extensions) # , raw_css=[css])

    pn.config.js_files.update({
        # This for copy to clipboard.
        "clipboard": "https://cdn.jsdelivr.net/npm/clipboard@2/dist/clipboard.min.js",
    })

    #pn.config.js_files.update({
    #    'ngl': 'https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js',
    #})
    #pn.extension(comms='ipywidgets')

    return pn


def gen_id(n=1, pre="uuid-"):
    """
    Generate ``n`` universally unique identifiers prepended with ``pre`` string.
    Return string if n == 1 or list of strings if n > 1
    """
    # The HTML4 spec says:
    # ID and NAME tokens must begin with a letter ([A-Za-z]) and may be followed by any number of letters,
    # digits ([0-9]), hyphens ("-"), underscores ("_"), colons (":"), and periods (".").
    import uuid
    if n == 1:
        return pre + str(uuid.uuid4())
    elif n > 1:
        return [pre + str(uuid.uuid4()) for i in range(n)]
    else:
        raise ValueError("n must be > 0 but got %s" % str(n))


def get_template_cls_from_name(name):
    """
    Return panel template from string.
    Support name in the form `FastList` as well as `FastListTemplate`.
    """
    # Example: pn.template.FastGridTemplate or pn.template.GoldenTemplate
    if hasattr(pn.template, name):
        return getattr(pn.template, name)

    try_name =  name + "Template"
    if hasattr(pn.template, try_name):
        return getattr(pn.template, try_name)

    raise ValueError(f"""
Don't know how to return panel template from string: {name}
Possible templates are: {list(pn.template.__dict__.keys())}
""")


class HTMLwithClipboardBtn(pn.pane.HTML):
    """
    Receives an HTML string and returns an HTML pane with a button that allows the user
    to copy the content to the system clipboard.
    Requires call to abipanel to load JS the extension.
    """

    # This counter is shared by all the instances. We use it so that the js script is included only once.
    _init_counter = [0]

    def __init__(self, object=None, btn_cls=None, **params):
        super().__init__(object=object, **params)

        self._init_counter[0] += 1
        my_id = gen_id()
        btn_cls = "bk bk-btn bk-btn-default" if btn_cls is None else str(btn_cls)

        # Build new HTML string with js section if first call.
        new_text = f"""
<div id="{my_id}">{self.object}</div>
<br>
<button class="clip-btn {btn_cls}" type="button" data-clipboard-target="#{my_id}"> Copy to clipboard </button>
<hr>
"""
        if self._init_counter[0] == 1:
            new_text += "<script> $(document).ready(function() {new ClipboardJS('.clip-btn')}) </script> "

        self.object = new_text


def mpl(fig, sizing_mode='stretch_width', with_controls=False, with_divider=True, **kwargs):
    """
    Helper function returning a panel Column with a matplotly pane followed by
    a divider and (optionally) controls to customize the figure.
    """
    #try:
    #    from ipympl.backend_nbagg import FigureManager, Canvas, is_interactive
    #    interactive = True
    #except:
    #    interactive = False

    col = pn.Column(sizing_mode=sizing_mode); ca = col.append
    mpl_pane = pn.pane.Matplotlib(fig, **kwargs)
    ca(mpl_pane)

    if with_controls:
        ca(pn.Accordion(("matplotlib controls", mpl_pane.controls(jslink=True))))
        ca(pn.layout.Divider())

    if with_divider:
        ca(pn.layout.Divider())

    return col


def ply(fig, sizing_mode='stretch_width', with_chart_studio=True, with_help=True,
        with_divider=True, with_controls=False):
    """
    Helper function returning a panel Column with a plotly pane,  buttons to push the figure
    to plotly chart studio and, optionally, controls to customize the figure.
    """
    col = pn.Column(sizing_mode=sizing_mode); ca = col.append
    plotly_pane = pn.pane.Plotly(fig, config={'responsive': True})
    ca(plotly_pane)

    if with_chart_studio:
        md = pn.pane.Markdown("""

The button on the left allows you to **upload the figure** to the
plotly [chart studio server](https://plotly.com/chart-studio-help/)
so that it is possible to share the figure or customize it via the chart studio editor.
In order to use this feature, you need to create a free account following the instructions
reported [in this page](https://plotly.com/chart-studio-help/how-to-sign-up-to-plotly/).
Then **add the following section** to your $HOME/.pmgrc.yaml configuration file:

```yaml
PLOTLY_USERNAME: john_doe  # Replace with your username
PLOTLY_API_KEY: secret     # To get your api_key go to: profile > settings > regenerate key
```

so that AbiPy can authenticate your user on the chart studio portal before pushing the figure to the cloud.
If everything is properly configured, a new window is automatically created in your browser.
""")

        btn = pnw.Button(name="Upload to chart studio server")
        def push_to_cs(event):
            push_to_chart_studio(fig)
        btn.on_click(push_to_cs)

        if with_help:
            acc = pn.Accordion(("What is this?", md))
            ca(pn.Row(btn, acc))
        else:
            ca(pn.Row(btn))

    if with_controls:
        ca(pn.Accordion(("plotly controls", plotly_pane.controls(jslink=True))))

    if with_divider:
        ca(pn.layout.Divider())

    return col


def dfc(df,
        wdg_type="dataframe",
        #wdg_type="tabulator",  # More recent version
        with_export_btn=True, with_controls=False, transpose=False, **kwargs):
    """
    Helper function returning a panel Column with a DataFrame or Tabulator widget followed by
    a divider and (optionally) controls to customize the figure.

    Note that not all the options work as exected. See comments below.
    """
    if "disabled" not in kwargs: kwargs["disabled"] = True
    if "sizing_mode" not in kwargs: kwargs["sizing_mode"] = "stretch_width"
    #if "sizing_mode" not in kwargs: kwargs["sizing_mode"] = "scale_width"
    if transpose:
        df = df.transpose()

    if wdg_type == "dataframe":
        if "auto_edit" not in kwargs: kwargs["auto_edit"] = False
        w = pnw.DataFrame(df, **kwargs)
    elif wdg_type == "tabulator":
        # This seems to be buggy
        w = pnw.Tabulator(df, **kwargs) #theme="modern",
    else:
        raise ValueError(f"Don't know how to handle widget type: `{wdg_type}`")

    col = pn.Column(sizing_mode='stretch_width'); ca = col.append
    ca(w)

    if with_export_btn:

        # Define callbacks with closure.
        #clip_button = pnw.Button(name="Copy to clipboard")
        #def to_clipboard(event):
        #    df.to_clipboard()
        #clip_button.on_click(to_clipboard)

        def to_xlsx():
            """
            Based on https://panel.holoviz.org/gallery/simple/file_download_examples.html
            NB: This requires xlsxwriter package else pandas raises ModuleNotFoundError.
            """
            output = io.BytesIO()
            writer = pd.ExcelWriter(output, engine='xlsxwriter')
            df.to_excel(writer, sheet_name="DataFrame")
            writer.save()  # Important!
            output.seek(0)
            return output

        def to_latex():
            """Convert DataFrame to latex string."""
            output = io.StringIO()
            df.to_latex(buf=output)
            output.seek(0)
            return output

        def to_md():
            """Convert DataFrame to markdown string."""
            output = io.StringIO()
            df.to_markdown(buf=output)
            output.seek(0)
            return output

        def to_json():
            """Convert DataFrame to json string."""
            output = io.StringIO()
            df.to_json(path_or_buf=output)
            output.seek(0)
            return output

        d = dict(
            xlsx=pnw.FileDownload(filename="data.xlsx", callback=to_xlsx),
            tex=pnw.FileDownload(filename="data.tex", callback=to_latex),
            md=pnw.FileDownload(filename="data.md", callback=to_md),
            json=pnw.FileDownload(filename="data.json", callback=to_json),
            )

        def download(event):
            print(f'Clicked menu item: "{event.new}"')
            file_download = d[event.new]
            print(file_download)
            #file_download._clicks = -1
            print("Calling transfer")
            file_download._transfer()
            #return file_download.callback()

        # FIXME: Menu button occpies less space but the upload does not work
        #menu_btn = pnw.MenuButton(name='Export to:', items=list(d.keys()))
        #menu_btn.on_click(download)
        #ca(menu_btn)

        # For the time being we use a Row with buttons.
        ca(pn.Row(*d.values(), sizing_mode="scale_width"))

    if with_controls:
        ca(pn.Accordion(("dataframe controls", w.controls(jslink=True))))

    ca(pn.layout.Divider())

    return col


class MyMarkdown(pn.pane.Markdown):
    """
    A Markdown pane renders the markdown markup language to HTML and
    displays it inside a bokeh Div model. It has no explicit
    priority since it cannot be easily be distinguished from a
    standard string, therefore it has to be invoked explicitly.
    """

    extensions = param.List(default=[
        # Extensions used by the superclass.
        "extra", "smarty", "codehilite",
        # My extensions
        #'pymdownx.arithmatex',
        #'pymdownx.details',
        #"pymdownx.tabbed",
    ],

        doc="""Markdown extension to apply when transforming markup."""
    )


class ButtonContext(object):
    """
    A context manager for buttons triggering computations on the server.

    This manager disables the button when we __enter__ and changes the name of the button to "running".
    It reverts to the initial state of the button one __exit__ is invoked, showing the Exception type
    in a "red" button if an exception is raised during the computation.

    This a very important tool because we need to disable the button when we start the computation
    to prevent the user from triggering multiple callbacks while the server is still working.
    At the same time, whathever happens in the callback, the button should go back to "clickable" mode
    when the callback returns so that the user can try to change the parameters and rerun.

    Note also that we want to provide some graphical feedback to the user if something goes wrong.
    At present we don't expose the python traceback on the client.
    It would be nice but we need panel machinery to do that.
    Moreover this is not the recommended approach for security reasons so we just change the "color"
    of the button and use the string representation of the exception as button name.
    """

    def __init__(self, btn):
        self.btn = btn
        self.prev_name, self.prev_type = btn.name, btn.button_type

    def __enter__(self):
        # Disable the button.
        self.btn.name = "Running ..."
        self.btn.button_type = "warning"
        self.btn.disabled = True
        return self.btn

    def __exit__(self, exc_type, exc_value, traceback):
        # First of all, reenable the button so that the user can stil interact with the GUI.
        self.btn.disabled = False

        if exc_type:
            # Exception --> signal to the user that something went wrong for 4 seconds.
            self.btn.name = str(exc_type)
            self.btn.button_type = "danger"
            import time
            time.sleep(3)

        # Back to the original button state.
        self.btn.name, self.btn.button_type = self.prev_name, self.prev_type

        # Don't handle the exception
        return None


class AbipyParameterized(param.Parameterized):
    """
    Base classs for AbiPy panels. Provide widgets for parameters supported by the differet AbiPy methods.
    and helper functions to perform typical operations when build dashboards from python.
    """

    verbose = param.Integer(0, bounds=(0, None), doc="Verbosity Level")
    mpi_procs = param.Integer(1, bounds=(1, None), doc="Number of MPI processes used for running Abinit")

    # mode = "webapp" This flag may be used to limit the user and/or decide the options that should
    # be exposed. For instance struct_viewer == "Vesta" does not make sense in webapp mode

    warning = pn.pane.Markdown(
"""
Please **refresh** the page using the refresh button of the browser if plotly figures are not shown.
""")

    @lazy_property
    def mpl_kwargs(self):
        """Default arguments passed to AbiPy matplotlib plot methods."""
        return dict(show=False, fig_close=True)

    def pws(self, *keys):
        """
        Helper function returning the list of parameters and widgets defined in self from a list of strings.
        Accepts also widget or parameter instances.
        """
        items, miss = [], []
        for k in keys:
            if isinstance(k, str):
                if k in self.param:
                    items.append(self.param[k])
                elif hasattr(self, k):
                    items.append(getattr(self, k))
                else:
                    miss.append(k)
            else:
                # Assume widget instance.
                items.append(k)

        if miss:
            raise ValueError("Cannot find `%s` in param or in attribute space" % str(miss))

        return items

    def helpc(self, method_name, extra_items=None):
        """
        Add accordion with a brief description and a warning after the button.
        The description of the tool is taken from the docstring of the callback.
        Return Column.
        """
        col = pn.Column(); ca = col.append

        acc = pn.Accordion(("Help", pn.pane.Markdown(getattr(self, method_name).__doc__)))

        if hasattr(self, "warning"):
            acc.append(("Warning", self.warning))

        if extra_items is not None:
            for name, attr in extra_items:
                acc.append((name, item))

        ca(pn.layout.Divider())
        ca(acc)

        return col

    @staticmethod
    def html_with_clipboard_btn(html_str, **kwargs):
        return HTMLwithClipboardBtn(html_str, **kwargs)

    @staticmethod
    def get_software_stack():
        """Return column with version of python packages in tabular format."""
        from abipy.abilab import software_stack
        return pn.Column("## Software stack:", dfc(software_stack(as_dataframe=True), with_export_btn=True),
                         sizing_mode="scale_width")

    def get_template_from_tabs(self, tabs, template):
        """
        This method receives a panel Tabs, include them in a template and return the panel template.
        """
        if template is None:
            return tabs

        cls = get_template_cls_from_name(template)

        kwargs = dict(
            # A title to show in the header. Also added to the document head meta settings and as the browser tab title.
            title=self.__class__.__name__,
            header_background="#ff8c00 ", # Dark orange
            #header_color="#ff8c00",
            #favicon (str): URI of favicon to add to the document head (if local file, favicon is base64 encoded as URI).
            #logo (str): URI of logo to add to the header (if local file, logo is base64 encoded as URI).
            #sidebar_footer (str): Can be used to insert additional HTML. For example a menu, some additional info, links etc.
            #enable_theme_toggle=False,  # If True a switch to toggle the Theme is shown. Default is True.
        )

        template = cls(**kwargs)
        if hasattr(template.main, "append"):
            template.main.append(tabs)
        else:
            # Assume .main area acts like a GridSpec
            template.main[:,:] = tabs

        return template


class HasStructureParams(AbipyParameterized):
    """
    Mixin class for panel objects providing a |Structure| object.
    """
    # Viewer widgets.
    struct_view_btn = pnw.Button(name="View structure", button_type='primary')
    struct_viewer = pnw.Select(name="Select viewer", value="vesta",
                               options=["jsmol", "vesta", "xcrysden", "vtk", "crystalk", "ngl",
                                        "matplotlib", "ase_atoms", "mayavi"])

    @property
    def structure(self):
        """Structure object provided by the subclass."""
        raise NotImplementedError(f"Subclass {type(self)} should implement `structure` attribute.")

    @pn.depends("struct_view_btn.clicks")
    def view_structure(self):
        """Visualize input structure."""
        if self.struct_view_btn.clicks == 0: return

        with ButtonContext(self.struct_view_btn):
            v = self.struct_viewer.value

            if v == "jsmol":
                #pn.extension(comms='ipywidgets') #, js_files=js_files)
                view = self.structure.get_jsmol_view()
                #from ipywidgets_bokeh import IPyWidget
                #view = IPyWidget(widget=view) #, width=800, height=300)
                #import ipywidgets as ipw
                #from IPython.display import display
                #display(view)
                #return pn.Row(display(view))
                #view = pn.ipywidget(view)
                #view = pn.panel(view)
                #view = pn.pane.IPyWidget(view)
                #print(view)
                #view = pn.Column(view, sizing_mode='stretch_width')
                return view

            if v == "crystalk":
                view = self.structure.get_crystaltk_view()
                return pn.panel(view)

            if v == "ngl":
                #js_files = {'ngl': 'https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js'}
                #pn.extension(comms='ipywidgets', js_files=js_files)
                view = self.structure.get_ngl_view()
                #return pn.panel(view)

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
            #return pn.Row(applet)

            if v == "ase_atoms":
                return mpl(self.structure.plot_atoms(rotations="default", **self.mpl_kwargs))

            return self.structure.visualize(appname=self.struct_viewer.value)

    def get_struct_view_tab_entry(self):
        """
        Return tab entry to visualize the structure.
         """
        return ("Structure", pn.Row(
            pn.Column("# Visualize structure",
                      *self.pws("struct_viewer", "struct_view_btn", self.helpc("view_structure"))),
                       pn.Column(self.get_structure_info(),
                                 self.view_structure)
        ))

    def get_structure_info(self):
        """
        Column with lattice parameters, angles and atomic positions grouped by type.
        """
        return get_structure_info(self.structure)


def get_structure_info(structure):
    col = pn.Column(sizing_mode='scale_width'); ca = col.append; cext = col.extend

    d = structure.get_dict4pandas(with_spglib=True)

    keys = index = [#"formula", "natom", "volume",
                    "abi_spg_number",
                    "spglib_symb", "spglib_num",  "spglib_lattice_type"]
    df_spg = pd.Series(data=d, index=index).to_frame()
    cext(["# Spacegroup:", dfc(df_spg, with_export_btn=False)])

    # Build dataframe with lattice lenghts.
    rows = []; keys = ("a", "b", "c")
    rows.append({k: d[k] * abu.Ang_Bohr for k in keys})
    rows.append({k: d[k] for k in keys})
    df_len = pd.DataFrame(rows, index=["Å", "Bohr"]).transpose().rename_axis("Lattice lenghts")

    # Build dataframe with lattice angles.
    rows = []; keys =  ("alpha", "beta", "gamma")
    rows.append({k: d[k] for k in keys})
    rows.append({k: np.radians(d[k]) for k in keys})
    df_ang = pd.DataFrame(rows, index=["Degrees", "Radians"]).transpose().rename_axis("Lattice angles")

    cext(["# Lattice lengths:", dfc(df_len, with_export_btn=False)])
    cext(["# Lattice angles:", dfc(df_ang, with_export_btn=False)])
    #row = pn.Row(dfc(df_len, with_export_btn=False), dfc(df_ang, with_export_btn=False), sizing_mode="scale_width")
    #ca(row)

    # Build dataframe with atomic positions grouped by element symbol.
    symb2df = structure.get_symb2coords_dataframe()
    accord = pn.Accordion(sizing_mode='stretch_width')
    for symb, df in symb2df.items():
        accord.append((f"Coordinates of {symb} sites:", dfc(df, with_export_btn=False)))
    ca(accord)

    return col


class PanelWithNcFile(AbipyParameterized):
    """
    This frame allows the user to inspect the dimensions and the variables reported in a netcdf file.
    Tab showing information on the netcdf file.

    Subclasses should implement the `ncfile` property
    """

    @property
    def ncfile(self):
        """abc does not play well with parametrized so we rely on this to enforce the protocol."""
        raise NotImplementedError("subclass should implement the `ncfile` property.")

    def get_ncfile_panel(self):
        col = pn.Column(sizing_mode='stretch_width'); ca = col.append

        #nc_grpname = pnw.Select(name="nc group name", options=["/"])

        # Get dataframe with dimensions.
        dims_df = self.ncfile.get_dims_dataframe(path="/")
        ca(dfc(dims_df))

        #vars_df = self.ncfile.get_dims_dataframe(path="/")
        #ca(dfc(vars_df))

        #ca(("NCFile", pn.Row(
        #    pn.Column("# NC dimensions and variables",
        #              dfc(dims_df, wdg_type="tabulator"),
        #             )),
        #))

        return col


class PanelWithElectronBands(AbipyParameterized):
    """
    Mixin class for panel object associated to AbiPy object providing an |ElectronBands| object.

    Subclasses should implement the `ebands` property.
    """

    # Bands plot
    with_gaps = pnw.Checkbox(name='show gaps')
    #ebands_ylims
    #ebands_e0
    # e0: Option used to define the zero of energy in the band structure plot. Possible values:
    #     - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
    #     -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
    #     -  None: Don't shift energies, equivalent to e0=0
    set_fermie_to_vbm = pnw.Checkbox(name="Set Fermie to VBM")

    plot_ebands_btn = pnw.Button(name="Plot e-bands", button_type='primary')

    # e-DOS plot.
    edos_method = pnw.Select(name="e-DOS method", options=["gaussian", "tetra"])
    edos_step = pnw.Spinner(name='e-DOS step (eV)', value=0.1, step=0.05, start=1e-6, end=None)
    edos_width = pnw.Spinner(name='e-DOS Gaussian broadening (eV)', value=0.2, step=0.05, start=1e-6, end=None)
    plot_edos_btn = pnw.Button(name="Plot e-DOS", button_type='primary')

    # SKW interpolation of KS band energies.
    skw_lpratio = pnw.IntInput(name='lpratio', value=5, step=1, start=1, end=None)
    skw_line_density = pnw.IntInput(name='skw_linedensity', value=20, step=1, start=1, end=None)
    plot_skw_btn = pnw.Button(name="Plot SKW", button_type='primary')

    # For the max size of file see: https://github.com/holoviz/panel/issues/1559
    ebands_kpath = None
    #ebands_kpath_fileinput = pnw.FileInput(accept=".nc")
    ebands_kpath_fileinput = param.FileSelector()
    ebands_kmesh = None
    #ebands_kmesh_fileinput = pnw.FileInput(accept=".nc")
    ebands_kmesh_fileinput = param.FileSelector()

    # Fermi surface plotter.
    #fs_viewer = pnw.Select(name="FS viewer", options=["matplotlib", "xcrysden"])
    #plot_fermi_surface_btn = pnw.Button(name="Plot Fermi surface", button_type='primary')

    @property
    def ebands(self):
        """abc does not play well with parametrized so we rely on this to enforce the protocol."""
        raise NotImplementedError("subclass should implement `ebands` property.")

    @pn.depends("ebands_kpath_fileinput", watch=True)
    def get_ebands_kpath(self):
        """
        """
        print(type(self.ebands_kpath_fileinput))
        #print(self.ebands_kpath_fileinput)
        #print(self.ebands_kpath_fileinput.path)
        #if self.ebands_kpath_fileinput.value is None: return None
        #filename = self.ebands_kpath_fileinput
        bdata = self.ebands_kpath_fileinput

        workdir = tempfile.mkdtemp()
        fd, tmp_path = tempfile.mkstemp(suffix=".nc") #filename)
        print(tmp_path)
        with open(tmp_path, "wb") as fh:
            fh.write(bdata)
            from abipy.electrons import ElectronBands
            self.ebands_kpath = ElectronBands.from_file(tmp_path)
            return self.ebands_kpath

    def get_plot_ebands_widgets(self):
        """Column with the widgets used to plot ebands."""
        return pn.Column(self.with_gaps, self.set_fermie_to_vbm, self.plot_ebands_btn)

    @pn.depends('plot_ebands_btn.clicks')
    def on_plot_ebands_btn(self):
        """Button triggering ebands plot."""
        if self.plot_ebands_btn.clicks == 0: return

        with ButtonContext(self.plot_ebands_btn):
            if self.set_fermie_to_vbm.value:
                self.ebands.set_fermie_to_vbm()

            #fig1 = self.ebands.plot(e0="fermie", ylims=None, with_gaps=self.with_gaps.value, max_phfreq=None,
            #                        fontsize=8, **self.mpl_kwargs)

            col = pn.Column(sizing_mode='stretch_width'); ca = col.append
            ca("## Electronic band structure:")
            fig1 = self.ebands.plotly(e0="fermie", ylims=None, with_gaps=self.with_gaps.value, max_phfreq=None,
                                      show=False)
            ca(ply(fig1))

            ca("## Brillouin zone and k-path:")
            #kpath_pane = mpl(self.ebands.kpoints.plot(**self.mpl_kwargs), with_divider=False)
            kpath_pane = ply(self.ebands.kpoints.plotly(show=False), with_divider=False)
            df_kpts = self.ebands.kpoints.get_highsym_datataframe()
            ca(pn.Row(kpath_pane, df_kpts))
            ca(pn.layout.Divider())

            #ca(bkw.PreText(text=self.ebands.to_string(verbose=self.verbose)))

            return col

    def get_plot_edos_widgets(self):
        """Widgets to compute the e-DOS."""
        return pn.Column(self.edos_method, self.edos_step, self.edos_width, self.plot_edos_btn)

    @pn.depends('plot_edos_btn.clicks')
    def on_plot_edos_btn(self):
        """Button triggering edos plot."""
        if self.plot_edos_btn.clicks == 0: return

        with ButtonContext(self.plot_edos_btn):
            edos = self.ebands.get_edos(method=self.edos_method.value,
                                        step=self.edos_step.value, width=self.edos_width.value)

            return pn.Row(mpl(edos.plot(**self.mpl_kwargs)), sizing_mode='scale_width')
            #return pn.Row(ply(edos.plotly(show=False))), sizing_mode='scale_width')

    @pn.depends('plot_skw_btn.clicks')
    def on_skw_btn(self):
        """
        """
        if self.plot_skw_btn.clicks == 0: return
        print("ebands_kpath:", self.ebands_kpath)
        #print("ebands_kpath.filename:", self.ebands_kpath_fileinput.filename)

        with ButtonContext(self.plot_skw_btn):
            col = pn.Column(sizing_mode='stretch_width'); ca = col.append

            d = self.ebands.interpolate(lpratio=self.skw_lpratio.value, line_density=self.skw_line_density.value,
                                        #vertices_names=None
                                        kmesh=None, is_shift=None, bstart=0, bstop=None, filter_params=None, verbose=0)

            ca("## SKW interpolated bands")
            ca(ply(d.ebands_kpath.plotly(with_gaps=self.with_gaps.value, show=False)))
            if self.ebands_kpath is not None:
                ca(ply(self.ebands_kpath.plotly(show=False)))

            return col

    def get_plot_skw_widgets(self):
        """Widgets to compute e-DOS."""

        foo = self.ebands_kmesh_fileinput = pn.Param(
            self.param['ebands_kpath_fileinput'],
            widgets={'ebands_kpath_fileinput': pn.widgets.FileInput}
        )

        return pn.Row(pn.Column(
            self.skw_lpratio, self.skw_line_density, self.with_gaps, foo,
            #self.ebands_kmesh_fileinput,
            self.plot_skw_btn),
            self.on_skw_btn)

    #def get_plot_fermi_surface_widgets(self):
    #    """Widgets to compute e-DOS."""
    #    return pn.Column(self.fs_viewer, self.plot_fermi_surface_btn)

    #@pn.depends('plot_fermi_surface_btn.clicks')
    #def on_plot_fermi_surface_btn(self):
    #    if self.plot_fermi_surface_btn.clicks == 0: return

    #    # Cache eb3d
    #    if hasattr(self, "_eb3d"):
    #        eb3d = self._eb3d
    #    else:
    #        # Build ebands in full BZ.
    #        eb3d = self._eb3d = self.ebands.get_ebands3d()

    #    if self.fs_viewer.value == "matplotlib":
    #        # Use matplotlib to plot isosurfaces corresponding to the Fermi level (default)
    #        # Warning: requires skimage package, rendering could be slow.
    #        fig = eb3d.plot_isosurfaces(e0="fermie", cmap=None, **self.mpl_kwargs)
    #        return pn.Row(mpl(fig), sizing_mode='scale_width')

    #    else:
    #        raise ValueError("Invalid choice: %s" % self.fs_viewer.value)

    #    #elif self.fs_viewer.value == "xcrysden":
    #        # Alternatively, it's possible to export the data in xcrysden format
    #        # and then use `xcrysden --bxsf mgb2.bxsf`
    #        #eb3d.to_bxsf("mgb2.bxsf")
    #        # If you have mayavi installed, try:
    #        #eb3d.mvplot_isosurfaces()



class BaseRobotPanel(AbipyParameterized):
    """Base class for panels with AbiPy robot."""

    compare_params_btn = pnw.Button(name="Compare structures", button_type='primary')
    transpose_params = pnw.Checkbox(name='Transpose tables')

    @pn.depends("compare_params_btn.clicks")
    def on_compare_params_btn(self):
        """
        """
        if self.compare_params_btn.clicks == 0: return

        with ButtonContext(self.compare_params_btn):
            col = pn.Column(sizing_mode='stretch_width'); ca = col.append
            transpose = self.transpose_params.value

            dfs = self.robot.get_structure_dataframes()
            ca("# Lattice daframe")
            ca(dfc(dfs.lattice, transpose=transpose))

            ca("# Params daframe")
            ca(dfc(self.robot.get_params_dataframe(), transpose=transpose))

            accord = pn.Accordion(sizing_mode='stretch_width')
            accord.append(("Atomic positions", dfc(dfs.coords, transpose=transpose)))
            ca(accord)

            return col

    # TODO: widgets to change robot labels.

    def get_compare_params_widgets(self):
        """
        """
        row = pn.Row(pn.Column(
            self.compare_params_btn, self.transpose_params),
            self.on_compare_params_btn,
            sizing_mode="scale_both")

        return row


class PanelWithEbandsRobot(BaseRobotPanel):
    """
    Mixin class for panels with a robot that owns a list of of |ElectronBands|
    """

    # Widgets to plot ebands.
    ebands_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot",
        options=["gridplot", "combiplot", "boxplot", "combiboxplot"]) # "animate",
    ebands_plotter_btn = pnw.Button(name="Plot", button_type='primary')
    ebands_df_checkbox = pnw.Checkbox(name='With Ebands DataFrame', value=False)

    # Widgets to plot edos.
    edos_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot",
        options=["gridplot", "combiplot"])
    edos_plotter_btn = pnw.Button(name="Plot", button_type='primary')

    def get_ebands_plotter_widgets(self):
        return pn.Column(self.ebands_plotter_mode, self.ebands_df_checkbox, self.ebands_plotter_btn)

    @pn.depends("ebands_plotter_btn.clicks")
    def on_ebands_plotter_btn(self):
        if self.ebands_plotter_btn.clicks == 0: return

        with ButtonContext(self.ebands_plotter_btn):

            ebands_plotter = self.robot.get_ebands_plotter()
            plot_mode = self.ebands_plotter_mode.value
            plot_func = getattr(ebands_plotter, plot_mode, None)
            if plot_func is None:
                raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)

            fig = plot_func(**self.mpl_kwargs)
            col = pn.Column(mpl(fig), sizing_mode='scale_width')
            if self.ebands_df_checkbox.value:
                df = ebands_plotter.get_ebands_frame(with_spglib=True)
                col.append(dfc(df))

            return pn.Row(col, sizing_mode='scale_width')

    def get_edos_plotter_widgets(self):
        return pn.Column(self.edos_plotter_mode, self.edos_plotter_btn)

    @pn.depends("edos_plotter_btn.clicks")
    def on_edos_plotter_btn(self):
        """Plot the electronic density of states."""
        if self.edos_plotter_btn.clicks == 0: return

        with ButtonContext(self.edos_plotter_btn):
            edos_plotter = self.robot.get_edos_plotter()
            plot_mode = self.edos_plotter_mode.value
            plot_func = getattr(edos_plotter, plot_mode, None)
            if plot_func is None:
                raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)

            fig = plot_func(**self.mpl_kwargs)

            return pn.Row(pn.Column(mpl(fig)), sizing_mode='scale_width')
