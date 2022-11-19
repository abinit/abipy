""""Basic tools and mixin classes for AbiPy panels."""
from __future__ import annotations

import io
import sys
import tempfile
import functools
import textwrap
import time
import traceback
import shutil
import param
import numpy as np
import pandas as pd
import bokeh.models.widgets as bkw
import panel as pn
import panel.widgets as pnw

from monty.functools import lazy_property
from monty.termcolor import cprint
from abipy.core import abinit_units as abu
from abipy.core.structure import Structure
from abipy.tools.plotting import push_to_chart_studio
from abipy.tools.decorators import Appender


_ABINIT_TEMPLATE_NAME = "FastList"


def set_abinit_template(template_name):
    global _ABINIT_TEMPLATE_NAME
    _ABINIT_TEMPLATE_NAME = template_name


def get_abinit_template_cls_kwds():
    cls =  get_template_cls_from_name(_ABINIT_TEMPLATE_NAME)
    kwds = dict(header_background="#ff8c00", # Dark orange
                favicon="/assets/img/abinit_favicon.ico",
                logo="/assets/img/abipy_logo.png", # TODO: Need new panel version to fix logo alignment in FastLIst.
                #sidebar_footer (str): Can be used to insert additional HTML.
                #                      For example a menu, some additional info, links etc.
                #enable_theme_toggle=False,  # If True a switch to toggle the Theme is shown. Default is True.
                )

    return cls, kwds


def open_html(html_string: str, browser: str = None):
    """
    Open a string with an HTML document in browser.
    """
    import tempfile
    import webbrowser
    with tempfile.NamedTemporaryFile(mode="wt", suffix=".html", delete=False) as tmpfile:
        tmpfile.write(html_string)
        webbrowser.get(browser).open_new_tab(f"file://{tmpfile.name}")


def abipanel(panel_template: str = "FastList"):
    """
    Activate panel extensions used by AbiPy. Return panel module.

    Args:
        panel_template: String with the name of the panel template to be used by default.
    """
    try:
        import panel as pn
    except ImportError as exc:
        cprint("Use `conda install panel` or `pip install panel` to install the python package.", "red")
        raise exc

    set_abinit_template(panel_template)

    #pn.extension(loading_spinner='dots', loading_color='#00aa41')

    pn.extension(notifications=True)
    #pn.config.notifications = True
    #pn.state.notifications.position = 'top-right'

    extensions = [
        "plotly",
        #"katex",
        "mathjax",
        "terminal",
        "tabulator",
        "ace",   # NB: This enters in conflict with Abipy Book
        #"gridstack",
        #"ipywidgets",
    ]

    css_files = [
        # FIXME
        #pn.io.resources.CSS_URLS['font-awesome'],
    ]

    #pn.extension(loading_spinner='petal', loading_color='#00aa41')
    #print("loading extensions:", extensions)

    #import os
    #abipy_css = os.path.join(os.path.dirname(__file__), "assets", "css", "abipy.css")

    pn.extension(*extensions, css_files=css_files) #, css_files=[abipy_css])
    #pn.extension(template='fast', theme='dark')

    pn.config.js_files.update({
        # This for copy to clipboard.
        "clipboard": "https://cdn.jsdelivr.net/npm/clipboard@2/dist/clipboard.min.js",
        # This for the jsmol viewer.
        "jsmol": "https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js",
    })

    #pn.extension('ipywidgets')

    #pn.config.js_files.update({
    #    'ngl': 'https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js',
    #})
    #pn.extension(comms='ipywidgets')
    #pn.config.sizing_mode = "stretch_width"

    #pn.config.css_files.extend([
    #    "https://maxcdn.bootstrapcdn.com/font-awesome/4.5.0/css/font-awesome.min.css",
    #])

    #css = """
    #    .pnx-file-upload-area input[type=file] {
    #        width: 100%;
    #        height: 100%;
    #        border: 3px dashed #9E9E9E;
    #        background: transparent;
    #        border-radius: 5px;
    #        text-align: left;
    #        margin: auto;
    #    }
    #"""

    #css = open(abipy_css, "rt").read()
    #pn.config.raw_css.append(css)

    return pn


def gen_id(n: int = 1, pre: str = "uuid-"):
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


def get_template_cls_from_name(name: str):
    """
    Return panel template from string.
    Support name in the form `FastList` as well as `FastListTemplate`.
    """
    # Example: pn.template.FastGridTemplate or pn.template.GoldenTemplate
    if hasattr(pn.template, name):
        return getattr(pn.template, name)

    try_name = name + "Template"
    if hasattr(pn.template, try_name):
        return getattr(pn.template, try_name)

    raise ValueError(f"""
Don't know how to return panel template from string: `{name}`.
Possible templates are: {list(pn.template.__dict__.keys())}
""")


add_mp_rest_docstring = Appender("""
NB: Requires a String API key for accessing the MaterialsProject REST interface.
Please apply on the Materials Project website for one.
If this is None, the code will check if there is a `PMG_MAPI_KEY` in your .pmgrc.yaml.
If so, it will use that environment
This makes easier for heavy users to simply add this environment variable to their setups and MPRester can
then be called without any arguments.
""", indents=0)


def depends_on_btn_click(btn_name: str,
                         show_doc: bool = True,
                         show_shared_wdg_warning: bool = True,
                         show_exc: bool = True):
    """
    This decorator is used for callbacks triggered by a button of name `btn_name`

    Args:
        btn_name: String with the name of the button.
        show_doc: If True, a Markdown pane with the doc string is returned the first time.
        show_shared_warning:
        show_exc: If True, a Markdown pane with the backtrace is returned if an exception is raised.
    """
    def decorator(func):
        @functools.wraps(func)
        def decorated(*args, **kwargs):
            self = args[0]
            btn = getattr(self, btn_name)
            if btn.clicks == 0:
                if show_doc and self._enable_show_doc:
                    doc = getattr(self, func.__name__).__doc__
                    if doc is None:
                        doc = f"\nNo docstring found for function `{func.__name__}`\n"
                    doc = textwrap.dedent(doc)
                    doc = f"## Description\n\n{doc}\n\n"
                    #print(doc)
                    objects = [my_md(doc)]
                    if show_shared_wdg_warning and self._enable_show_shared_wdg_warning:
                        warning = pn.pane.Alert(SHARED_WIDGETS_WARNING, alert_type="danger")
                        objects.append(warning)
                    return pn.Column(*objects, sizing_mode="stretch_width")
                else:
                    return None

            with ButtonContext(btn):
                return func(*args, **kwargs)

        f = pn.depends(f"{btn_name}.clicks")(decorated)
        if show_exc: f = show_exception(f)
        return f

    return decorator


def show_exception(func):
    """
    This decorator returns a Markdown pane with the backtrace
    if the function raises an exception.
    """
    @functools.wraps(func)
    def decorated(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as exc:
            s = traceback.format_exc()
            print(s)
            return pn.pane.Markdown(f"```shell\n{s}\n```", sizing_mode="stretch_both")

    return decorated


class HTMLwithClipboardBtn(pn.pane.HTML):
    """
    Receives an HTML string and returns an HTML pane with a button that allows users
    to copy the content to the system clipboard.
    Requires call to abipanel to load the JS extension.
    """

    # This counter is shared by all the instances.
    # We use it so that the js script is included only once.
    #_init_counter = [0]

    def __init__(self, object=None, btn_cls=None, **params):
        super().__init__(object=object, **params)

        #self._init_counter[0] += 1
        my_id = gen_id()
        btn_cls = "bk bk-btn bk-btn-default" if btn_cls is None else str(btn_cls)

        # Build new HTML string with js section if first call.
        new_text = f"""
<div id="{my_id}">{self.object}</div>
<br>
<button class="clip-btn {btn_cls}" type="button" data-clipboard-target="#{my_id}"> Copy to clipboard </button>
<hr>

"""
        if True: # self._init_counter[0] == 1:
            #new_text += " <script> $(document).ready(function() {new ClipboardJS('.clip-btn')}) </script> "
            # $(document).ready(function() {
            new_text += """ <script>
if (typeof abipy_clipboard === 'undefined') {
    var abipy_clipboard = new ClipboardJS('.clip-btn');
}

if (typeof abipy_notyf === 'undefined') {
    // Create an instance of Notyf
    var abipy_notyf = new Notyf();

    abipy_clipboard.on('success', function(e) {
        abipy_notyf.success('Text copied to clipboard');
    });

    abipy_clipboard.on('error', function(e) {
        abipy_notyf.error('Cannot copy text to clipboard');
    });
}
</script> """

        self.object = new_text


def mpl(fig, sizing_mode='stretch_width', with_controls=False, with_divider=True, **kwargs) -> pn.Column:
    """
    Helper function returning a panel Column with a matplotly pane followed by
    a divider and (optionally) controls to customize the figure.
    """
    col = pn.Column(sizing_mode=sizing_mode); ca = col.append

    #try:
    #    import ipympl
    #    has_ipympl = True
    #except ImportError:
    #    has_ipympl = False

    #if "interactive" not in kwargs and has_ipympl: # "ipympl" in sys.modules:
    #    print("mpl in interactive mode")
    #    kwargs["interactive"] = True

    if "tight" not in kwargs: kwargs["tight"] = True

    mpl_pane = pn.pane.Matplotlib(fig, **kwargs)
    ca(mpl_pane)

    if with_controls:
        ca(pn.Accordion(("matplotlib controls", mpl_pane.controls(jslink=True))))
        ca(pn.layout.Divider())

    if with_divider:
        ca(pn.layout.Divider())

    return col


def ply(fig, sizing_mode='stretch_both', with_chart_studio=False, with_help=False,
        with_divider=True, with_controls=False) -> pn.Column:
    """
    Helper function returning a panel Column with a plotly pane,  buttons to push the figure
    to plotly chart studio and, optionally, controls to customize the figure.
    """
    col = pn.Column(sizing_mode=sizing_mode); ca = col.append

    config = dict(
      responsive=True,
      #showEditInChartStudio=True,
      showLink=True,
      plotlyServerURL="https://chart-studio.plotly.com",
      )

    plotly_pane = pn.pane.Plotly(fig, config=config)
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

        btn = pnw.Button(name="Upload figure to chart studio server")

        def push_to_cs(event):
            with ButtonContext(btn):
                push_to_chart_studio(fig)

        btn.on_click(push_to_cs)

        if with_help:
            acc = pn.Accordion(("Help", md))
            ca(pn.Row(btn, acc))
        else:
            ca(pn.Row(btn))

    if with_controls:
        ca(pn.Accordion(("plotly controls", plotly_pane.controls(jslink=True))))

    if with_divider:
        ca(pn.layout.Divider())

    return col


def dfc(df: pd.DataFrame,
        #wdg_type: str = "dataframe",
        wdg_type: str ="tabulator",  # More recent version. Still problematic
        with_export_btn=True, with_controls=False, with_divider=True, transpose=False, **kwargs):
    """
    Helper function returning a panel Column with a DataFrame or Tabulator widget followed by
    a divider and (optionally) controls to customize the figure.

    Note that not all the options work as exected. See comments below.
    """
    if "disabled" not in kwargs: kwargs["disabled"] = True
    #if "sizing_mode" not in kwargs: kwargs["sizing_mode"] = "stretch_both"
    if "sizing_mode" not in kwargs: kwargs["sizing_mode"] = "scale_width"
    if transpose:
        df = df.transpose()

    if wdg_type == "dataframe":
        if "auto_edit" not in kwargs: kwargs["auto_edit"] = False
        w = pnw.DataFrame(df, **kwargs)
    elif wdg_type == "tabulator":
        w = pnw.Tabulator(df, **kwargs)
    else:
        raise ValueError(f"Don't know how to handle widget type: `{wdg_type}`")

    col = pn.Column(sizing_mode='stretch_both'); ca = col.append
    ca(w)

    if with_export_btn:
        # Define callbacks with closure.
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

        # For the time being we use a Row with buttons.
        #ca(pn.Row(*d.values(), sizing_mode="scale_width"))
        ca(pn.Card(*d.values(), title="Export table", collapsed=True,
                   sizing_mode='stretch_width', header_color="blue",
        ))

        #def download(event):
        #    file_download = d[event.new]
        #    #print(f'Clicked menu item: "{event.new}"')
        #    print(file_download)
        #    #file_download._clicks = -1
        #    file_download._transfer()
        #    return file_download.callback()

        # FIXME: Menu button occupies less space but the upload does not work
        #menu_btn = pnw.MenuButton(name='Export table', items=list(d.keys()))
        #menu_btn.on_click(download)
        #ca(menu_btn)

    if with_controls:
        ca(pn.Accordion(("dataframe controls", w.controls(jslink=True))))

    if with_divider:
        ca(pn.layout.Divider())

    return col


def my_md(string: str, **kwargs) -> pn.pane.Markdown:
    """
    Return a Markdown pane from `string`.
    Extra kwargs are passed to pane.Markdown.
    The string can contain links to Abinit variables in the wikilink format.

    .. example::

        The anaddb variable [[dipdip@anaddb]] has the same meaning as the Abinit variable [[dipdip]].
    """
    import re
    WIKILINK_RE = r'\[\[([^\[]+)\]\]'
    from abipy.abio.abivar_database.variables import get_codevars
    vars_code = get_codevars()

    def repl(match):
        var_name = match.group(1).strip()
        codename = "abinit"
        i = var_name.find("@")
        if i != -1:
            var_name, codename = var_name[:i], var_name[i+1:]
        try:
            var = vars_code[codename][var_name]
            return var.html_link()
        except:
            return f"WRONG LINK for {var_name}"

    string = re.sub(WIKILINK_RE, repl, string)

    # TODO: Latex
    extra_extensions = [
         "markdown.extensions.admonition",
         #'pymdownx.arithmatex',
         #'pymdownx.details',
         #"pymdownx.tabbed",
    ]

    md = pn.pane.Markdown(string, **kwargs)
    md.extensions.extend(extra_extensions)

    return md


class ButtonContext():
    """
    A context manager for buttons triggering computations on the server.

    This manager disables the button when we __enter__ and changes the name of the button to "running".
    It reverts to the initial state of the button once __exit__ is invoked, showing the Exception type
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

    def __init__(self, btn: pnw.Button):
        self.btn = btn
        self.prev_name, self.prev_type = btn.name, btn.button_type

    def __enter__(self):
        # Disable the button.
        self.btn.name = "Running ..."
        self.btn.button_type = "warning"
        self.btn.loading = True
        self.btn.disabled = True
        return self.btn

    def __exit__(self, exc_type, exc_value, traceback):
        # First of all, reenable the button so that the user can stil interact with the GUI.
        self.btn.disabled = False
        self.btn.loading = False

        if exc_type:
            # Exception --> signal to the user that something went wrong for 2 seconds.
            self.btn.name = str(exc_type)
            self.btn.button_type = "danger"
            time.sleep(2)

        # Back to the original button state.
        self.btn.name, self.btn.button_type = self.prev_name, self.prev_type


class Loading():
    """
    A context manager for setting the loading attribute of a panel object.
    """

    def __init__(self, panel_obj, err_wdg=None, width=70):
        self.panel_obj = panel_obj
        self.err_wdg = err_wdg
        if err_wdg is not None: self.err_wdg.object = ""
        self.width = int(width)

    def __enter__(self):
        self.panel_obj.loading = True
        return self.panel_obj

    def __exit__(self, exc_type, exc_value, traceback):
        self.panel_obj.loading = False

        if self.err_wdg is not None:
            if exc_type:
                self.err_wdg.object = "```sh\n%s\n```" % textwrap.fill(str(exc_value), width=self.width)


class ActiveBar():
    """
    A context manager that sets progress.active to True on entry and False when we exit.
    """

    def __init__(self, progress, err_wdg=None, width=70):
        self.progress = progress
        self.err_wdg = err_wdg
        if err_wdg is not None: self.err_wdg.object = ""
        self.width = int(width)

    def __enter__(self):
        self.progress.active = True
        return self.progress

    def __exit__(self, exc_type, exc_value, traceback):
        self.progress.active = False

        if exc_type:
            # Change the color to signal the user that something bad happened.
            self.progress.bar_color = "danger"

            if self.err_wdg is not None:
                from textwrap import fill
                self.err_wdg.object = "```sh\n%s\n```" % fill(str(exc_value), width=self.width)


SHARED_WIDGETS_WARNING = """
Note some of the widgets are **shared by the different tabs**.
If you change the value of a variables in the active tab,
the same value will **automagically** appear in all the other tabs sharing the same
widget yet results/figures are not automatically recomputed."""

#In other words, if you change some variable in the active tab and then you move to another tab,
#the results/figures (if any) are stil computed with the **old input** hence you will have to
#recompute the new results by clicking the button."""


class AbipyParameterized(param.Parameterized):
    """
    Base class for AbiPy panels. Provides helper functions for typical operations needed for
    building dashboard and basic parameters supported by the subclasses.
    """

    verbose = param.Integer(0, bounds=(0, None), doc="Verbosity Level")
    mpi_procs = param.Integer(1, bounds=(1, None), doc="Number of MPI processes used for running Fortran code.")

    plotly_template = param.ObjectSelector(default="plotly",
                                           objects=["plotly", "plotly_white", "plotly_dark", "ggplot2",
                                                    "seaborn", "simple_white", "none"])

    # This flag is set to True if we are serving apps from the Abinit server.
    # It is used to impose limitations on what users can do and select the options that should be exposed.
    # For instance, structure_viewer == "Vesta" does not make sense in we are not serving from a local server.
    #
    has_remote_server = param.Boolean(False)

    warning = pn.pane.Markdown(SHARED_WIDGETS_WARNING, name="warning")

    def __init__(self, **params):

        super().__init__(**params)
        if self.has_remote_server:
            self.param.mpi_procs.bounds = (1, 1)
            print("Changing mpi_procs.bounds")
            print("self.param.mpi_procs:", self.param.mpi_procs.bounds)

        # This internal variables are used to enable/disable the output of the docstring and the warning message.
        self._enable_show_doc = True
        self._enable_show_shared_wdg_warning = True

    #def abiopen(self, filepath):
    #    from abipy.abilab import abiopen
    #    try:
    #        abiobj = abiopen(filepath)
    #        return abiobj
    #    finally:
    #        self._abiopened_files.append(filepath)

    def disable_show_messages(self, show_doc=False, show_shared_wdg_warning=False):
        """
        Context manager to disable the output of the docstring and the warning message.

        with self.disable_show_messages():
            pn.Row("## Callback with no info message", self.on_plot_phbands_and_phdos)

        """
        from abipy.tools.context_managers import temporary_change_attributes
        return temporary_change_attributes(self,
                                           _enable_show_doc=show_doc,
                                           _enable_show_shared_wdg_warning=show_shared_wdg_warning,
                                           )

    @pn.depends("plotly_template")
    def on_plotly_template(self):
        """
        Change the default plotly template.
        NB: It's not a good idea to expose this option when running the server as
        other users will be affected by this change hence this function is just for
        internal use.
        """
        import plotly.io as pio
        pio.templates.default = self.plotly_template

    @lazy_property
    def mpl_kwargs(self) -> dict:
        """Default arguments passed to AbiPy matplotlib plot methods."""
        return dict(show=False, fig_close=True)

    def pws_col(self, keys, **kwargs) -> pn.Column:
        return pn.Column(*self.pws(keys), **kwargs)

    def pws_row(self, keys, **kwargs) -> pn.Row:
        return pn.Row(*self.pws(keys), **kwargs)

    def wdg_box(self, keys, **kwargs) -> pn.WidgetBox:
        return pn.WidgetBox(*self.pws(keys), **kwargs)

    def pws(self, keys):
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
                elif k.startswith("#"):
                    # Markdown string
                    items.append(k)
                    #items.append(pn.layout.Divider(margin=(-20, 0, 0, 0)))
                else:
                    miss.append(k)
            else:
                # Assume widget instance.
                items.append(k)

        #for item in items: print("item", item, "of type:", type(item))
        if miss:
            raise ValueError(f"Cannot find `{str(miss)}` in param or in attribute space")

        return items

    def get_summary_view_for_abiobj(self, abiobj, **kwargs):
        text = abiobj.to_string(verbose=self.verbose)

        view = pnw.Terminal(output=f"\n\n{text}",
            #height=1200, # Need this one else the terminal is not show properly
            sizing_mode='stretch_both',
        )
        #view = pn.Row(bkw.PreText(text=text, sizing_mode="scale_both"))
        return view

    def wdg_exts_with_get_panel(self, name='File extensions supported:'):
        """
        Return Select widget with the list of file extensions implementing a get_panel method.
        """
        from abipy.abilab import extcls_supporting_panel
        exts = [e[0] for e in extcls_supporting_panel(as_table=False)]
        return pnw.Select(name=name, options=exts)

    @staticmethod
    def html_with_clipboard_btn(html_str: str, **kwargs):
        if hasattr(html_str, "_repr_html_"):
            html_str = html_str._repr_html_()

        return HTMLwithClipboardBtn(html_str, **kwargs)

    @staticmethod
    def get_software_stack() -> pn.Column:
        """Return column with version of python packages in tabular format."""
        from abipy.abilab import software_stack
        return pn.Column("## Software stack:",
                         #pn.layout.Divider(margin=(-20, 0, 0, 0)),
                         dfc(software_stack(as_dataframe=True), with_export_btn=False),
                         sizing_mode="scale_width",
                         )

    @staticmethod
    def get_fileinput_section(file_input) -> pn.Column:
        # All credits go to:
        # https://github.com/MarcSkovMadsen/awesome-panel/blob/master/application/pages/styling/fileinput_area.py
        #
        css_style = """
        <style>
        .pnx-file-upload-area input[type=file] {
            width: 100%;
            height: 100%;
            border: 3px dashed #9E9E9E;
            background: transparent;
            border-radius: 5px;
            text-align: left;
            margin: auto;
        }
        </style>"""

        return pn.Column(pn.pane.HTML(css_style, width=0, height=0, sizing_mode="stretch_width", margin=0),
                         file_input, sizing_mode="stretch_width")

    @staticmethod
    def get_abifile_from_file_input(file_input, use_structure=False):
        #print("filename", file_input.filename, "\nvalue", file_input.value)
        workdir = tempfile.mkdtemp()

        fd, tmp_path = tempfile.mkstemp(suffix=file_input.filename)
        with open(tmp_path, "wb") as fh:
            fh.write(file_input.value)

        from abipy.abilab import abiopen
        abifile = abiopen(tmp_path)
        if use_structure:
            abifile = Structure.as_structure(abifile)
            # Remove the file since it's not needed anymore.
            shutil.rmtree(tmp_path, ignore_errors=True)

        return abifile

    def get_ebands_from_file_input(self, file_input, remove=True):
        """
        Read and return an |ElectronBands| object from a file_input widget.
        Return None if the file does not provide an ebands object.
        Remove the file if remove == True.
        """
        with self.get_abifile_from_file_input(file_input) as abifile:
            ebands = getattr(abifile, "ebands", None)
        if remove: abifile.remove()
        return ebands

    @staticmethod
    def get_alert_data_transfer() -> pn.pane.Alert:
        # https://discourse.holoviz.org/t/max-upload-size/2121/5
        return pn.pane.Alert("""
Please note that this web interface is not designed to handle **large data transfer**.
To post-process the data stored in a big file e.g. a WFK.nc file,
we strongly suggest executing the **abigui.py**  script on the same machine where the file is hosted.
Also, note that examples of post-processing scripts are available in the
[AbiPy gallery](https://abinit.github.io/abipy/gallery/index.html).

Last but not least, keep in mind that **the file extension matters** when uploading a file
so don't change the default extension used by ABINIT.
Also, use `.abi` for ABINIT input files and `.abo` for the main output file.
""", alert_type="info")

    @staticmethod
    def get_template_cls_from_name(template):
        return get_template_cls_from_name(template)

    def get_abinit_template_cls_kwds(self):
        return get_abinit_template_cls_kwds()

    def get_template_from_tabs(self, tabs, template, **tabs_kwargs):
        """
        This method receives panel Tabs or a dictionary,
        include them in a template and return the template instance.
        """
        if isinstance(tabs, dict):
            if "sizing_mode" not in tabs_kwargs: tabs_kwargs["sizing_mode"] = "stretch_width"
            tabs = pn.Tabs(*tabs.items(), **tabs_kwargs)

        if template is None or str(template) == "None":
            #tabs.append(return_to_top_html())
            return tabs

        cls = get_template_cls_from_name(template)
        #cls, kwargs = get_abinit_template_cls_kwds()

        kwargs = dict(
            # A title to show in the header. Also added to the document head meta settings and as the browser tab title.
            title=self.__class__.__name__,
            header_background="#ff8c00", # Dark orange
            favicon="/assets/img/abinit_favicon.ico",
            logo="/assets/img/abipy_logo.png", # TODO: Need new panel version to fix logo alignment in FastLIst.
            #sidebar_footer (str): Can be used to insert additional HTML. For example a menu, some additional info, links etc.
            #enable_theme_toggle=False,  # If True a switch to toggle the Theme is shown. Default is True.
        )

        template = cls(**kwargs)
        if hasattr(template.main, "append"):
            template.main.append(tabs)
            #template.main.append(return_to_top_html())
        else:
            # Assume main area acts like a GridSpec
            template.main[:,:] = tabs

        return template


def return_to_top_html():
    print("In return to top")

    #<a id="button"></a>

    html = r"""
<!--
Back to top button
https://codepen.io/matthewcain/pen/ZepbeR
-->

<style>
#button {
  display: inline-block;
  background-color: #FF9800;
  width: 50px;
  height: 50px;
  text-align: center;
  border-radius: 4px;
  position: fixed;
  bottom: 30px;
  right: 30px;
  transition: background-color .3s,
    opacity .5s, visibility .5s;
  opacity: 0;
  visibility: hidden;
  z-index: 1000;
}
#button::after {
  content: "\f077";
  font-family: FontAwesome;
  font-weight: normal;
  font-style: normal;
  font-size: 2em;
  line-height: 50px;
  color: #fff;
}
#button:hover {
  cursor: pointer;
  background-color: #333;
}
#button:active {
  background-color: #555;
}
#button.show {
  opacity: 1;
  visibility: visible;
}
</style>

<script>

// Code executed on page ready
$(function() {

    var btn = $('#button');
    console.log("button", btn);

    $(window).scroll(function() {
      if ($(window).scrollTop() > 300) {
        console.log("show");
        btn.addClass('btn show');
      } else {
        btn.removeClass('show');
        console.log("btn hide");
      }
    });

    btn.on('click', function(e) {
      e.preventDefault();
      $('html, body').animate({scrollTop:0}, '300');
    });

})
</script>
"""

    return pn.pane.HTML(html, width=0, height=0)


class PanelWithStructure(AbipyParameterized):
    """
    A paremeterized object with a |Structure| object.
    """

    structure_viewer = param.ObjectSelector(default="jsmol",
                                            objects=["jsmol", "vesta", "xcrysden", "vtk", "crystalk", "ngl",
                                                     "matplotlib", "plotly", "ase_atoms", "mayavi"])

    def __init__(self, structure: Structure, **params):

        super().__init__(**params)
        self.structure = structure

        if self.has_remote_server:
            # Change the list of allowed visualizers.
            self.param.structure_viewer.objects = ["jsmol", "crystalk", "ngl", "matplotlib", "plotly", "ase_atoms"]

        self.view_structure_btn = pnw.Button(name="View structure", button_type='primary')

    @depends_on_btn_click('view_structure_btn', show_shared_wdg_warning=False)
    def on_view_structure(self):
        """Visualize input structure."""
        v = self.structure_viewer

        if v == "jsmol":
            return jsmol_html(self.structure)

            #pn.extension(comms='ipywidgets') #, js_files=js_files)
            #view = self.structure.get_jsmol_view()
            #from ipywidgets_bokeh import IPyWidget
            #view = IPyWidget(widget=view) #, width=800, height=300)
            #import ipywidgets as ipw
            #from IPython.display import display
            #display(view)
            #return pn.Row(display(view))
            #return pn.ipywidget(view)
            #return pn.panel(view)
            #return pn.pane.IPyWidget(view)
            #print(view)
            #view = pn.Column(view, sizing_mode='stretch_width')
            #return view

        if v == "crystalk":
            view = self.structure.get_crystaltk_view()
            return pn.panel(view)

        if v == "plotly":
            return ply(self.structure.plotly(show=False))

        if v == "ngl":
            from pymatgen.io.babel import BabelMolAdaptor
            from pymatgen.io.xyz import XYZ
            # string_data = self.structure.to(fmt="xyz")

            #writer = BabelMolAdaptor(self)
            #string_data = str(XYZ(self.structure))
            #adapt = BabelMolAdaptor.from_string(string_data, file_format="xyz")
            ##pdb_string =
            #print(pdb_string)

            #from awesome_panel_extesions.pane.widgets.ngl_viewer import NGLViewer
            #view = NGLViewer()

            view.pdb_string = pdb_string
            return view

            #js_files = {'ngl': 'https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js'}
            #pn.extension(comms='ipywidgets', js_files=js_files)
            #view = self.structure.get_ngl_view()
            #return pn.panel(view)

            #pn.config.js_files["ngl"]="https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js"
            #pn.extension()

            html = """<div id="viewport" style="width:100%; height:100%;"></div>
            <script>
            stage = new NGL.Stage("viewport");
            stage.loadFile("rcsb://1NKT.mmtf", {defaultRepresentation: true});
            </script>"""

            ngl_pane = pn.pane.HTML(html, height=500, width=500)
            return pn.Row(ngl_pane)
            view = self.structure.get_ngl_view()

        #return self.structure.crystaltoolkitview()
        #import nglview as nv
        #view = nv.demo(gui=False)

        if v == "ase_atoms":
            return mpl(self.structure.plot_atoms(rotations="default", **self.mpl_kwargs))

        return self.structure.visualize(appname=self.structure_viewer)

    def get_structure_view(self) -> pn.Row:
        """
        Return Row with widgets to visualize the structure.
        """
        return pn.Row(
            self.pws_col(["## Visualize structure",
                          "structure_viewer",
                          "view_structure_btn",
                          ]),
            pn.Column(self.on_view_structure, self.get_structure_info())
        )

    def get_structure_info(self) -> pn.Column:
        """
        Return Column with lattice parameters, angles and atomic positions grouped by type.
        """
        return get_structure_info(self.structure)


def get_structure_info(structure: Structure) -> pn.Column:
    """
    Return Column with lattice parameters, angles and atomic positions grouped by type.
    """
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
    rows = []; keys = ("alpha", "beta", "gamma")
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


class NcFileViewer(AbipyParameterized):
    """
    This class implements toool to inspect dimensions and variables stored in a netcdf file.

    Relyes on the API provided by `AbinitNcFile` defined in core.mixins.py
    """

    nc_path = param.String("/", doc="nc group")

    def __init__(self, ncfile, **params):
        super().__init__(**params)
        self.ncfile = ncfile
        self.netcdf_info_btn = pnw.Button(name="Show info", button_type='primary')

    def get_ncfile_view(self) -> pn.Column:
        return pn.Column(
                self.netcdf_info_btn,
                self.on_netcdf_info_btn,
                sizing_mode='stretch_width',
        )

    @depends_on_btn_click('netcdf_info_btn')
    def on_netcdf_info_btn(self) -> pn.Column:
        """
        This Tab allows one to
        """
        # TODO: Finalize the implementation.
        col = pn.Column(sizing_mode='stretch_width'); ca = col.append

        #nc_grpname = pnw.Select(name="nc group name", options=["/"])
        input_string = self.ncfile.get_input_string()
        ca(f"## Input String")
        ca(bkw.PreText(text=input_string))

        #ca(f"## Global attributes")

        # Get dataframe with dimensions.
        dims_df = self.ncfile.get_dims_dataframe(path=self.nc_path)
        ca(f"## Dimensions in nc group: {self.nc_path}")
        ca(dfc(dims_df))
        #ca(f"## Variables")

        return col


class PanelWithElectronBands(PanelWithStructure):
    """
    Provide widgets and views for operating on |ElectronBands| object.
    """

    # Bands plot
    with_gaps = param.Boolean(False)
    with_kpoints_plot = param.Boolean(False)

    #ebands_ylims
    #ebands_e0
    # e0: Option used to define the zero of energy in the band structure plot. Possible values:
    #     - `fermie`: shift all eigenvalues to have zero energy at the Fermi energy (`self.fermie`).
    #     -  Number e.g e0=0.5: shift all eigenvalues to have zero energy at 0.5 eV
    #     -  None: Don't shift energies, equivalent to e0=0
    set_fermie_to_vbm = param.Boolean(False, label="Set Fermi energy to VBM")

    # e-DOS plot.
    edos_method = param.ObjectSelector(default="gaussian", label="Integration method for e-DOS",
                                       objects=["gaussian", "tetra"])
    edos_step_ev = param.Number(0.1, bounds=(1e-6, None), step=0.1, label='e-DOS step in eV')
    edos_width_ev = param.Number(0.2, step=0.05, bounds=(1e-6, None), label='e-DOS Gaussian broadening in eV')

    # SKW interpolation of the KS band energies (Abipy version).
    skw_lpratio = param.Integer(5, bounds=(1, None), label="SKW lpratio")
    skw_line_density = param.Integer(20, label="SKW line density")
    skw_ebands_kpath = None
    skw_ebands_kpath_fileinput = param.FileSelector(path="*.nc")

    # Fermi surface.
    fs_viewer = param.ObjectSelector(default="matplotlib", objects=["matplotlib", "xcrysden"])

    # Ifermi UI.
    ifermi_wigner_seitz = param.Boolean(True, label="Use Wigner Seitz cell",
                                        doc="Controls whether the cell is the Wigner-Seitz cell" +
                                            "or the reciprocal unit cell parallelepiped.")
    ifermi_interpolation_factor = param.Integer(default=8, label="Interpolation factor", bounds=(1, None),
                                                doc="The factor by which the band structure will be interpolated.")

    ifermi_eref = param.ObjectSelector(default="fermie", label="Energy reference",
                                       objects=["fermie", "cbm", "vbm"])

    ifermi_with_velocities = param.Boolean(True, label="Compute group velocities",
            doc="Generate the Fermi surface and calculate the group velocity at the center of each triangular face")
    ifermi_offset_eV = param.Number(default=0.0, label="Energy offset (eV) from energy reference",
                                    doc="Energy offset from the Fermi energy at which the isosurface is calculated.")
    ifermi_plot_type = param.ObjectSelector(default="plotly", label="Plot type", objects=["plotly", "matplotlib"])

    # These are used to implement plots in which we need to upload an additional file
    # For instance bands + edos.
    # For the max size of file see: https://github.com/holoviz/panel/issues/1559
    ebands_kpath = None
    ebands_kpath_fileinput = param.FileSelector(path="*.nc")
    ebands_kmesh = None
    ebands_kmesh_fileinput = param.FileSelector(path=".nc")

    effmass_accuracy = param.Integer(default=4, bounds=(1, None), label="Finite difference accuracy")
    effmass_degtol_ev = param.Number(default=1e-3, bounds=(0, None), label="Window in eV above/below the CBM/VBM")
    effmass_spin = param.ObjectSelector(default=0, objects=[0, 1], label="Spin index")

    def __init__(self, ebands, **params):

        self.ebands = ebands
        PanelWithStructure.__init__(self, structure=ebands.structure, **params)

        # Create buttons
        self.plot_ebands_btn = pnw.Button(name="Plot e-bands", button_type='primary')
        self.plot_edos_btn = pnw.Button(name="Plot e-DOS", button_type='primary')
        self.plot_skw_btn = pnw.Button(name="Plot SKW interpolant", button_type='primary')

        # Fermi surface plotter.
        #objects = [None, "matplotlib", "xcrysden"]
        #if self.has_remote_server: objects = [None, "matplotlib"]
        #self.plot_fs_viewer_btn = pnw.Button(name="Plot SKW interpolant", button_type='primary')

        self.plot_ifermi_btn = pnw.Button(name="Plot Fermi surface", button_type='primary')
        #self.ifermi_plane_normal = pnw.LiteralInput(name='Plane normal (list)', value=[0, 0, 0], type=list,
        #                                            placeholder="Enter normal in reduced coordinates")
        #self.ifermi_distance = pn.widgets.RangeSlider(
        #        name='distance', start=0, end=2 * max(ebands.structure.reciprocal_lattice.abc),
        #        value=(0, 0), step=0.01)

        #ebands_kpath_fileinput = pnw.FileInput(accept=".nc")
        #ebands_kmesh_fileinput = pnw.FileInput(accept=".nc")

        self.plot_effmass_btn = pnw.Button(name="Plot effective masses", button_type='primary')
        if ebands.nsppol != 2:
            self.param.effmass_spin.objects = [0]

    @staticmethod
    def _get_ebands_from_bstring(bstring):
        from abipy.electrons import ElectronBands
        return ElectronBands.from_binary_string(bstring)

    @pn.depends("ebands_kpath_fileinput", watch=True)
    def get_ebands_kpath(self):
        """
        Receives the netcdf file selected by the user as binary string.
        """
        self.ebands_kpath = self._get_ebands_from_bstring(self.ebands_kpath_fileinput)

    @pn.depends("ebands_kmesh_fileinput", watch=True)
    def get_ebands_kmesh(self):
        """
        Receives the netcdf file selected by the user as binary string.
        """
        self.ebands_kmesh = self._get_ebands_from_bstring(self.ebands_kmesh_fileinput)

    @pn.depends("skw_ebands_kpath_fileinput", watch=True)
    def get_skw_ebands_kpath(self):
        """
        Receives the netcdf file selected by the user as binary string.
        """
        self.skw_ebands_kpath = self._get_ebands_from_bstring(self.skw_ebands_kpath_fileinput)

    def get_plot_ebands_view(self) -> pn.Row:
        return pn.Row(
            self.pws_col(["### e-Bands Plot Options",
                          "with_gaps", "set_fermie_to_vbm", "with_kpoints_plot", "plot_ebands_btn",
                         ]),
            self.on_plot_ebands_btn
        )

    @depends_on_btn_click('plot_ebands_btn')
    def on_plot_ebands_btn(self) -> pn.Column:
        """
        This Tab allows one to plot the KS energies stored in the netcdf file
        as well as the associated list of **k**-points in the Brillouin.
        """
        if self.set_fermie_to_vbm:
            self.ebands.set_fermie_to_vbm()

        sz_mode = "stretch_width"
        col = pn.Column(sizing_mode=sz_mode); ca = col.append
        ca("## Electronic band structure:")
        fig1 = self.ebands.plotly(e0="fermie", ylims=None, with_gaps=self.with_gaps, max_phfreq=None, show=False)
        ca(ply(fig1))

        nkpt = len(self.ebands.kpoints)
        ktype = "IBZ sampling" if self.ebands.kpoints.is_ibz else "**k**-path"
        max_nkpt = 2000

        if self.with_kpoints_plot:
            ca(f"## Brillouin zone and {ktype}:")
            if nkpt < max_nkpt:
                kpath_pane = ply(self.ebands.kpoints.plotly(show=False), with_divider=False)
                df_kpts = dfc(self.ebands.kpoints.get_highsym_datataframe(), with_divider=False)
                ca(pn.Row(kpath_pane, df_kpts, sizing_mode=sz_mode))
            else:
                ca(f"k-points won't be shown as nkpt: {nkpt} is greater than {max_nkpt}")
        ca(pn.layout.Divider())

        return col

    def get_plot_edos_view(self) -> pn.Row:
        return pn.Row(
                self.pws_col(["## E-DOS Options", "edos_method", "edos_step_ev",
                              "edos_width_ev", "plot_edos_btn"]),
                self.on_plot_edos_btn
                )

    @depends_on_btn_click('plot_edos_btn')
    def on_plot_edos_btn(self) -> pn.Row:
        """
        Button triggering edos plot.
        """
        edos = self.ebands.get_edos(method=self.edos_method, step=self.edos_step_ev, width=self.edos_width_ev)

        return pn.Row(ply(edos.plotly(show=False)), sizing_mode='scale_width')

    def get_skw_view(self) -> pn.Row:
        """
        Column with widgets to use SKW.
        """
        wdg = pn.Param(
            self.param['skw_ebands_kpath_fileinput'],
            widgets={'skw_ebands_kpath_fileinput': pn.widgets.FileInput}
        )

        return pn.Row(
            self.pws_col(["## SKW options",
                          "skw_lpratio", "skw_line_density", "with_gaps",
                          "## Upload GSR.nc file with *ab-initio* energies along a k-path to compare with", wdg,
                          "plot_skw_btn",
                          ]),
            self.on_plot_skw_btn)

    @depends_on_btn_click('plot_skw_btn')
    def on_plot_skw_btn(self) -> pn.Column:
        """
        Button triggering SKW plot.
        """
        col = pn.Column(sizing_mode='stretch_width'); ca = col.append

        intp = self.ebands.interpolate(lpratio=self.skw_lpratio, line_density=self.skw_line_density,
                                       kmesh=None, is_shift=None, bstart=0, bstop=None, filter_params=None,
                                       verbose=self.verbose)

        ca("## SKW interpolated bands along an automatically selected high-symmetry **k**-path")
        ca(ply(intp.ebands_kpath.plotly(with_gaps=self.with_gaps, show=False)))

        if self.skw_ebands_kpath is not None:
            ca("## Input bands taken from file uploaded by user:")
            ca(ply(self.skw_ebands_kpath.plotly(with_gaps=self.with_gaps, show=False)))

            # Use line_density 0 to interpolate on the same set of k-points given in self.skw_ebands_kpath
            vertices_names = []
            for kpt in self.skw_ebands_kpath.kpoints:
                vertices_names.append((kpt.frac_coords, kpt.name))

            intp = self.ebands.interpolate(lpratio=self.skw_lpratio, vertices_names=vertices_names, line_density=0,
                                            kmesh=None, is_shift=None, bstart=0, bstop=None, filter_params=None,
                                            verbose=self.verbose)

            plotter = self.skw_ebands_kpath.get_plotter_with("Input", "Interpolated", intp.ebands_kpath)
            ca("## Input bands vs SKW interpolated bands:")
            ca(ply(plotter.combiplotly(show=False)))

        return col

    def get_effmass_view(self) -> pn.Row:
        """
        Return Row with widgets to compute effective masses with finite diff.
        """
        return pn.Row(
            self.pws_col(["### Effective masses options",
                          "effmass_accuracy",
                          "effmass_degtol_ev",
                          "effmass_spin",
                          "plot_effmass_btn",
                         ]),
            self.on_plot_effmass_btn
        )

    @depends_on_btn_click('plot_effmass_btn')
    def on_plot_effmass_btn(self) -> pn.Column:
        """
        Compute and visualize effective masses with finite differences.

        Note that the quality of the results strongly depend on the step i.e.
        the separation between two consecutive points along the k-path.
        The accuracy option allows one to change the number of points for the finite difference.
        """
        emana = self.ebands.get_effmass_analyzer()
        acc = self.effmass_accuracy
        degtol_ev = self.effmass_degtol_ev
        spin = self.effmass_spin

        col = pn.Column(sizing_mode='stretch_width'); ca = col.append

        if emana.select_vbm():
            ca(f"## Effective masses at the VBM with accuracy {acc}:")
            for segment in emana.segments:
                ca(mpl(segment.plot_emass(acc=acc, spin=spin, degtol_ev=degtol_ev, show=False)))
                ca("### effmass wrt accuracy and step: %.3f Ang-1" % segment.dk)
                df = segment.get_dataframe_with_accuracies()
                ca(dfc(df, with_export_btn=False))

        if emana.select_cbm():
            ca(f"## Effective masses at the CBM with accuracy {acc}:")
            for segment in emana.segments:
                ca(mpl(segment.plot_emass(acc=acc, spin=spin, degtol_ev=degtol_ev, show=False)))
                ca("### effmass wrt accuracy and step: %.3f Ang-1" % segment.dk)
                df = segment.get_dataframe_with_accuracies()
                ca(dfc(df, with_export_btn=False))

        return col

    def get_ifermi_view(self):
        """
        Widgets to visualize the Fermi surface with ifermi package.
        """
        controls = self.pws_col([
                          "ifermi_offset_eV", "ifermi_eref", "ifermi_wigner_seitz", "ifermi_interpolation_factor",
                          "ifermi_with_velocities", # "ifermi_plane_normal", "ifermi_distance",
                          "ifermi_plot_type", "plot_ifermi_btn",
                          ])

        return pn.Row(pn.Column("## ifermi options", controls), self.on_plot_ifermi_btn)

    @depends_on_btn_click('plot_ifermi_btn')
    def on_plot_ifermi_btn(self):
        """
        This Tab allows you to interpolate KS energies defined in the IBZ
        and visualize isosurfaces in the BZ using the [ifermi](https://fermisurfaces.github.io/IFermi/) package.

        The energy of the isosurface is given by: `energy_reference` + `energy_offset`.

        * use `fermie` with zero offset to visualize the Fermi surface in metals.
        * for electron pockets in semiconductors, use `cbm` with a positive offset e.g. 0.1 eV.
        * for hole pockets in semiconductors, use `vbm` and a  negative offset e.g. -0.1 eV.
        """
        # interpolate the energies onto a dense k-point mesh
        r = self.ebands.get_ifermi_fs(interpolation_factor=self.ifermi_interpolation_factor,
                                      mu=self.ifermi_offset_eV,
                                      eref=self.ifermi_eref,
                                      wigner_seitz=self.ifermi_wigner_seitz,
                                      calculate_dimensionality=False,
                                      with_velocities=self.ifermi_with_velocities,
                                      )

        plt = r.fs_plotter.get_plot(plot_type=self.ifermi_plot_type)

        if self.ifermi_plot_type == "plotly":
            fig = ply(plt, sizing_mode="stretch_both")
        elif self.ifermi_plot_type == "matplotlib":
            fig = plt.gcf()
            plt.close(fig=fig)
            fig = mpl(fig)
        else:
            raise NotImplementedError(f"{self.ifermi_plot_type}")

        col = pn.Column(sizing_mode="stretch_width")
        ca = col.append
        if self.ifermi_wigner_seitz:
            ca("## Energy isosurface in the Wigner-Seitz unit cell")
        else:
            ca("## Energy isosurface in the reciprocal unit cell parallelepiped")
        ca(fig)

        ene_range = [-3, 3]
        fig = self.ebands.boxplotly(e0=r.abs_isoenergy, ene_range=ene_range, show=False)
        ca(f"""
## Energy boxplot:

- Energy zero set at the absolute isoenergy: {r.abs_isoenergy:.3f} (eV)
- Energy range around zero: {ene_range} (eV)
""")
        ca(ply(fig, sizing_mode="stretch_both"))

        # TODO: This requires more testing
        #ifermi_plane_normal = self.ifermi_plane_normal.value
        #if any(c != 0 for c in ifermi_plane_normal):
        #    print("fs_plane normal:", ifermi_plane_normal)
        #    print("abc reciprocal", self.ebands.structure.reciprocal_lattice.abc)
        #    from ifermi.plot import FermiSlicePlotter
        #    for distance in self.ifermi_distance.value:
        #        # generate Fermi slice along the (0 0 1) plane going through the Γ-point.
        #        ca(f"## Fermi slice along the {ifermi_plane_normal} plane going through the Γ-point at distance: {distance}")
        #        fermi_slice = r.fs.get_fermi_slice(plane_normal=ifermi_plane_normal, distance=distance)
        #        slice_plotter = FermiSlicePlotter(fermi_slice)
        #        plt = slice_plotter.get_plot()
        #        fig = plt.gcf()
        #        plt.close(fig=fig)
        #        fig = mpl(fig)
        #        ca(fig)

        #ca(pn.layout.Divider())
        #ca("## Powered by [ifermi](https://fermisurfaces.github.io/IFermi/)")

        return col

    #@depends_on_btn_click('plot_fs_viewer_btn')
    #def on_plot_fs_viewer_btn(self):

    #    # Get eb3d (memoized)
    #    eb3d = self._eb3d = self.ebands.get_ebands3d()

    #    if self.fs_viewer == "matplotlib":
    #        # Use matplotlib to plot isosurfaces corresponding to the Fermi level (default)
    #        # Warning: requires skimage package, rendering could be slow.
    #        fig = eb3d.plot_isosurfaces(e0="fermie", cmap=None, **self.mpl_kwargs)
    #        return pn.Row(mpl(fig), sizing_mode='scale_width')

    #    elif self.fs_viewer == "xcrysden":
    #       # Alternatively, it's possible to export the data in xcrysden format
    #       # and then use `xcrysden --bxsf mgb2.bxsf`
    #       #eb3d.to_bxsf("mgb2.bxsf")
    #       # If you have mayavi installed, try:
    #       #eb3d.mvplot_isosurfaces()
    #       raise NotImplementedError(f"Invalid choice for fs_viewer: {self.fs_viewer}")

    #    else:
    #        raise ValueError(f"Invalid choice for fs_viewer: {self.fs_viewer}")

    #def get_fsviewer_view(self):
    #    """
    #    Widgets to visualize the Fermi surface with ifermi
    #    """
    #    #return pn.Row(pn.Column("## FS Viewer options", self.fs_viewer, plot_fs_viewer_btn),
    #    #                        self.on_plot_fs_viewer_btn)

    #    return pn.Row(
    #        self.pws_col(["## FS Viewer options",
    #                      "fs_viewer", "plot_fs_viewer_btn"
    #                      ]),
    #        self.on_plot_fs_viewer_btn)


class BaseRobotPanel(AbipyParameterized):
    """
    Base class for panels with AbiPy robot.
    """

    def __init__(self, robot, **params):
        self.robot = robot
        self.compare_params_btn = pnw.Button(name="Compare structures", button_type='primary')
        self.transpose_params = pnw.Checkbox(name='Transpose table', default=True)

        super().__init__(**params)

    @depends_on_btn_click("compare_params_btn")
    def on_compare_params_btn(self):
        """
        Compare lattice parameters and atomic positions.
        """
        col = pn.Column(sizing_mode='stretch_width'); ca = col.append
        transpose = self.transpose_params.value

        dfs = self.robot.get_structure_dataframes()
        ca("# Lattice dataframe")
        ca(dfc(dfs.lattice, transpose=transpose))

        ca("# Parameters dataframe")
        ca(dfc(self.robot.get_params_dataframe(), transpose=transpose))

        accord = pn.Accordion(sizing_mode='stretch_width')
        accord.append(("Atomic positions", dfc(dfs.coords, transpose=transpose)))
        ca(accord)

        return col

    # TODO: widgets to change robot labels.
    def get_compare_params_widgets(self):
        """
        """
        return pn.Row(pn.Column(
            self.compare_params_btn, self.transpose_params),
            self.on_compare_params_btn,
            sizing_mode="scale_both")


class PanelWithEbandsRobot(BaseRobotPanel):
    """
    Mixin class for panels with a robot that owns a list of of |ElectronBands|.
    """

    def __init__(self, robot, **params):

        BaseRobotPanel.__init__(self, robot=robot, **params)

        # Widgets to plot ebands.
        self.ebands_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot",
                                              options=["gridplot", "combiplot", "boxplot", "combiboxplot"]) # "animate",
        self.ebands_plotter_btn = pnw.Button(name="Plot", button_type='primary')
        self.ebands_df_checkbox = pnw.Checkbox(name='With Ebands DataFrame', value=False)

        # Widgets to plot edos.
        self.edos_plotter_mode = pnw.Select(name="Plot Mode", value="gridplot", options=["gridplot", "combiplot"])
        self.edos_plotter_btn = pnw.Button(name="Plot", button_type='primary')

    def get_ebands_plotter_widgets(self):
        return pn.Column(self.ebands_plotter_mode, self.ebands_df_checkbox, self.ebands_plotter_btn)

    @depends_on_btn_click("ebands_plotter_btn")
    def on_ebands_plotter_btn(self):
        """
        Plot the electronic density of states.
        """
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

    @depends_on_btn_click("edos_plotter_btn")
    def on_edos_plotter_btn(self):
        """
        Plot the electronic density of states.
        """
        edos_plotter = self.robot.get_edos_plotter()
        plot_mode = self.edos_plotter_mode.value
        plot_func = getattr(edos_plotter, plot_mode, None)
        if plot_func is None:
            raise ValueError("Don't know how to handle plot_mode: %s" % plot_mode)

        fig = plot_func(**self.mpl_kwargs)

        return pn.Row(pn.Column(mpl(fig)), sizing_mode='scale_width')


def jsmol_html(structure, supercell=(1, 1, 1), width=700, height=700, color="black", spin="false")  -> pn.Column:

    cif_str = structure.write_cif_with_spglib_symms(None, ret_string=True) #, symprec=symprec

    # There's a bug in boken when we use strings with several '" quotation marks
    # To bypass the problem I create a json list of strings and then I use a js variable
    # to recreate the cif file by joining the tokens with: var string = elements.join('\n');
    import json
    lines = cif_str.split("\n")
    lines = json.dumps(lines, indent=4)
    #print("lines:", lines)

    jsmol_div_id = gen_id()
    jsmol_app_name = "js1"
    supercell = "{" + ",".join(str(s) for s in supercell) + "}"
    #supercell = "{2, 2, 2}"

    # http://wiki.jmol.org/index.php/Jmol_JavaScript_Object/Functions#getAppletHtml

    html = f"""
<script type="text/javascript">
    $(document).ready(function() {{

       const lines = {lines};
       var script_str = 'load inline " ' + lines.join('\\n') + '" {supercell}';
       //console.log(script_str);

       var Info = {{
           color: "{color}",
           spin: {spin},
           antialiasDisplay: true,
           width: {width},
           height: {height},
           j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s",
           serverURL: "https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php",
           script: script_str,
           use: 'html5',
           disableInitialConsole: true,
           disableJ2SLoadMonitor: true,
           debug: false
       }};

       $("#{jsmol_div_id}").html(Jmol.getAppletHtml("{jsmol_app_name}", Info));
    }});
</script>

<div id="{jsmol_div_id}" style="height: 100%; width: 100%; position: relative;"></div>
"""

    #print(html)
    return pn.Column(pn.pane.HTML(html, sizing_mode="stretch_width"), sizing_mode="stretch_width")
