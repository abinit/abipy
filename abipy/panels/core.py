""""Basic tools and mixin classes for AbiPy panels."""

import param
import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.tools.plotting import push_to_chart_studio


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


class HTMLwithClipboardBtn(pn.pane.HTML):
    """
    Receives an HTML string and return an HTML pane with a button that allows the user
    to copy the content to the system clipboard.
    Requires call to abipanel to load JS extension.
    """

    # This counter is shared by all the instances. We use it so that the js script is included only once.
    _init_counter = [0]

    def __init__(self, object=None, **params):
        super().__init__(object=object, **params)

        self._init_counter[0] += 1
        my_id = gen_id()
        btn_cls = "bk bk-btn bk-btn-primary"

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


def mpl(fig, sizing_mode='stretch_width', with_controls=False, **kwargs):
    """
    Helper function returning a panel Column with a matplotly pane followed by
    a divider and (optionally) controls to customize the figure.
    """
    col = pn.Column(sizing_mode=sizing_mode); ca = col.append
    mpl_pane = pn.pane.Matplotlib(fig, **kwargs)
    ca(mpl_pane)
    ca(pn.layout.Divider())

    if with_controls:
        ca(pn.Accordion(("matplotlib controls", mpl_pane.controls(jslink=True))))
        ca(pn.layout.Divider())

    return col


def ply(fig, sizing_mode='stretch_width', with_chart_studio=True, with_controls=False):
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

        acc = pn.Accordion(("What is this?", md))

        button = pnw.Button(name="Upload to chart studio server", button_type='primary')
        def push_to_cs(event):
            push_to_chart_studio(fig)
        button.on_click(push_to_cs)

        ca(pn.Row(button, acc))

    if with_controls:
        ca(pn.Accordion(("plotly controls", plotly_pane.controls(jslink=True))))
        ca(pn.layout.Divider())

    ca(pn.layout.Divider())

    return col


def dfc(df, wdg_type="dataframe", with_copy_to_clipboard=True, with_controls=False, **kwargs):
    """
    Helper function returning a panel Column with a plotly pane followed by
    a divider and (optionally) controls to customize the figure.

    Note that not all the options work as exected. See comments below.
    """
    if "disabled" not in kwargs: kwargs["disabled"] = True
    if "sizing_mode" not in kwargs: kwargs["sizing_mode"] = "stretch_width"

    if wdg_type == "dataframe":
        if "auto_edit" not in kwargs: kwargs["auto_edit"] = False
        w = pnw.DataFrame(df, **kwargs)
    elif wdg_type == "tabulator":
        # This seems to be buggy
        w = pnw.Tabulator(df, **kwargs)
    else:
        raise ValueError(f"Don't know how to handle widget type: {wdg_type}")

    col = pn.Column(sizing_mode='stretch_width'); ca = col.append
    ca(w)

    if with_copy_to_clipboard:
        md = pn.pane.Markdown(r"""
The button on the left allows you to write a text representation of the dataframe to the system clipboard.
This can be pasted into Excel, for example.
"""
)
        acc = pn.Accordion(("What is this?", md))
        button = pnw.Button(name="Copy to clipboard", button_type='primary')
        def to_clipboard(event):
            df.to_clipboard()
        button.on_click(to_clipboard)

        ca(pn.Row(button, acc))
        #ca(pn.Row(button, md))

    if with_controls:
        # This seems to be buggy.
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
        'pymdownx.arithmatex',
        'pymdownx.details',
        "pymdownx.tabbed",
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

        # Back the original button state.
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


class HasStructureParams(AbipyParameterized):
    """
    Mixin class for panel objects providing a |Structure| object.
    """
    # Viewer widgets.
    struct_view_btn = pnw.Button(name="View structure", button_type='primary')
    struct_viewer = pnw.Select(name="Viewer", value="vesta",
                               options=["jsmol", "vesta", "xcrysden", "vtk", "crystalk", "ngl",
                                        "matplotlib", "ase_atoms", "mayavi"])

    @property
    def structure(self):
        """Structure object provided by the subclass."""
        raise NotImplementedError(f"Subclass {type(self)} should implement `structure` attribute.")

    def get_struct_view_tab_entry(self):
        """
        Return tab entry to visualize the structure.
         """
        return ("View Structure", pn.Row(
            pn.Column("# Visualize structure", *self.pws("struct_viewer", "struct_view_btn", self.helpc("view_structure"))),
            self.view_structure)
        )

    @param.depends("struct_view_btn.clicks")
    def view_structure(self):
        """Visualize input structure."""
        if self.struct_view_btn.clicks == 0: return

        with ButtonContext(self.struct_view_btn):
            v = self.struct_viewer.value

            if v == "jsmol":
                pn.extension(comms='ipywidgets') #, js_files=js_files)
                view = self.structure.get_jsmol_view()
                from ipywidgets_bokeh import IPyWidget
                view = IPyWidget(widget=view) #, width=800, height=300)
                #import ipywidgets as ipw
                from IPython.display import display
                #display(view)
                return pn.panel(view)
                #return pn.Row(display(view))

            if v == "crystalk":
                view = self.structure.get_crystaltk_view()
                return pn.panel(view)

            if v == "ngl":
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
            #return pn.Row(applet)

            if v == "ase_atoms":
                return mpl(self.structure.plot_atoms(rotations="default", **self.mpl_kwargs))

            return self.structure.visualize(appname=self.struct_viewer.value)


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

        # Get dataframe with dimesions.
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

    Subclasses should implement `ebands` property
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

    # Fermi surface plotter.
    fs_viewer = pnw.Select(name="FS viewer", options=["matplotlib", "xcrysden"])
    plot_fermi_surface_btn = pnw.Button(name="Plot Fermi surface", button_type='primary')

    @property
    def ebands(self):
        """abc does not play well with parametrized so we rely on this to enforce the protocol."""
        raise NotImplementedError("subclass should implement `ebands` property.")

    def get_plot_ebands_widgets(self):
        """Column with the widgets used to plot ebands."""
        return pn.Column(self.with_gaps, self.set_fermie_to_vbm, self.plot_ebands_btn)

    @param.depends('plot_ebands_btn.clicks')
    def on_plot_ebands_btn(self):
        """Button triggering ebands plot."""
        if self.plot_ebands_btn.clicks == 0: return

        with ButtonContext(self.plot_ebands_btn):
            if self.set_fermie_to_vbm.value:
                self.ebands.set_fermie_to_vbm()

            fig1 = self.ebands.plot(e0="fermie", ylims=None, with_gaps=self.with_gaps.value, max_phfreq=None,
                                    fontsize=8, **self.mpl_kwargs)

            fig2 = self.ebands.kpoints.plot(**self.mpl_kwargs)
            row = pn.Row(mpl(fig1), mpl(fig2)) #, sizing_mode='scale_width')
            text = bkw.PreText(text=self.ebands.to_string(verbose=self.verbose))

            return pn.Column(row, text, sizing_mode='scale_width')

    def get_plot_edos_widgets(self):
        """Widgets to compute e-DOS."""
        return pn.Column(self.edos_method, self.edos_step, self.edos_width, self.plot_edos_btn)

    @param.depends('plot_edos_btn.clicks')
    def on_plot_edos_btn(self):
        """Button triggering edos plot."""
        if self.plot_edos_btn.clicks == 0: return

        with ButtonContext(self.plot_edos_btn):
            edos = self.ebands.get_edos(method=self.edos_method.value,
                                        step=self.edos_step.value, width=self.edos_width.value)
            fig = edos.plot(**self.mpl_kwargs)
            return pn.Row(mpl(fig), sizing_mode='scale_width')

    def get_plot_fermi_surface_widgets(self):
        """Widgets to compute e-DOS."""
        return pn.Column(self.fs_viewer, self.plot_fermi_surface_btn)

    @param.depends('plot_fermi_surface_btn.clicks')
    def on_plot_fermi_surface_btn(self):
        if self.plot_fermi_surface_btn.clicks == 0: return

        # Cache eb3d
        if hasattr(self, "_eb3d"):
            eb3d = self._eb3d
        else:
            # Build ebands in full BZ.
            eb3d = self._eb3d = self.ebands.get_ebands3d()

        if self.fs_viewer.value == "matplotlib":
            # Use matplotlib to plot isosurfaces corresponding to the Fermi level (default)
            # Warning: requires skimage package, rendering could be slow.
            fig = eb3d.plot_isosurfaces(e0="fermie", cmap=None, **self.mpl_kwargs)
            return pn.Row(mpl(fig), sizing_mode='scale_width')

        else:
            raise ValueError("Invalid choice: %s" % self.fs_viewer.value)

        #elif self.fs_viewer.value == "xcrysden":
            # Alternatively, it's possible to export the data in xcrysden format
            # and then use `xcrysden --bxsf mgb2.bxsf`
            #eb3d.to_bxsf("mgb2.bxsf")
            # If you have mayavi installed, try:
            #eb3d.mvplot_isosurfaces()


class BaseRobotPanel(AbipyParameterized):
    """pass"""


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

    @param.depends("ebands_plotter_btn.clicks")
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

    @param.depends("edos_plotter_btn.clicks")
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
