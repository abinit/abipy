#!/usr/bin/env python

import base64
import datetime
import tempfile
import io
import os
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash_dangerously_set_inner_html import DangerouslySetInnerHTML

from dash.dependencies import Input, Output, State
from abipy import abilab
import abipy.data as abidata
from abipy.data.ucells import structure_from_ucell


def copy_to_clipboard_button(id):
    return DangerouslySetInnerHTML(f"""
<button class="btn btn-primary btn-sm" data-clipboard-target="#{id}"> Copy to clipboard </button>""")


def div_with_exc(exc):
    return html.Div([dbc.Jumbotron([html.H2("There was an error: %s" % str(exc),
                     className="text-danger")])])


app = dash.Dash(__name__,
                external_stylesheets=[dbc.themes.BOOTSTRAP],
                # external JavaScript files
                external_scripts=[
                    "https://cdn.jsdelivr.net/npm/clipboard@2/dist/clipboard.min.js",
                    #"https://code.jquery.com/jquery-3.4.1.min.js"
                ],
                # these meta_tags ensure content is scaled correctly on different devices
                # see: https://www.w3schools.com/css/css_rwd_viewport.asp for more
                meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}]
)

# Since we're adding callbacks to elements that don't exist in the app.layout,
# Dash will raise an exception to warn us that we might be doing something wrong.
# In this case, we're adding the elements through a callback, so we can ignore the exception.
app.config.suppress_callback_exceptions = True

# we use the Row and Col components to construct the sidebar header
# it consists of a title, and a toggle, the latter is hidden on large screens
sidebar_header = dbc.Row(
    [
        dbc.Col(html.H1("Actions", className="display-4")),
        dbc.Col(
            [
                html.Button(
                    # use the Bootstrap navbar-toggler classes to style
                    html.Span(className="navbar-toggler-icon"),
                    className="navbar-toggler",
                    # the navbar-toggler classes don't set color
                    style={
                        "color": "rgba(0,0,0,.5)",
                        "border-color": "rgba(0,0,0,.1)",
                    },
                    id="navbar-toggle",
                ),
                html.Button(
                    # use the Bootstrap navbar-toggler classes to style
                    html.Span(className="navbar-toggler-icon"),
                    className="navbar-toggler",
                    # the navbar-toggler classes don't set color
                    style={
                        "color": "rgba(0,0,0,.5)",
                        "border-color": "rgba(0,0,0,.1)",
                    },
                    id="sidebar-toggle",
                ),
            ],
            # the column containing the toggle will be only as wide as the
            # toggle, resulting in the toggle being right aligned
            width="auto",
            # vertically align the toggle in the center
            align="center",
        ),
    ]
)

pages = [
    dbc.NavLink("Upload Structure",  href="/page-1", id="page-1-link"),
    dbc.NavLink("Convert Structure", href="/page-2", id="page-2-link"),
    dbc.NavLink("Symmetry Analysis", href="/page-3", id="page-3-link"),
    dbc.NavLink("Kpath",             href="/page-4", id="page-4-link"),
    dbc.NavLink("Build GS Input",    href="/page-5", id="page-5-link"),
    #dbc.NavLink("Build Ebands Input",href="/page-5", id="page-5-link"),
    #dbc.NavLink("Build DFPT Input",href="/page-5", id="page-5-link"),
    #dbc.NavLink("Build G0W0 Input",href="/page-5", id="page-5-link"),
    dbc.NavLink("About",             href="/page-6", id="page-6-link"),
]
num_pages = len(pages)

sidebar = html.Div([
        sidebar_header,
        # we wrap the horizontal rule and short blurb in a div that can be hidden on a small screen
        html.Div([
            html.Hr(),
            #html.P("Actions:", className="lead"),
        ], id="blurb"),
        # use the Collapse component to animate hiding / revealing links
        dbc.Collapse(dbc.Nav(pages, vertical=True, pills=True), id="collapse"),
    ],
    id="sidebar",
)

#example_structure = structure_from_ucell("Si")

app.layout = html.Div([
    dcc.Location(id="url"),
    sidebar,
    html.Div(id="page-content"),
    # Hidden div inside the app that stores the intermediate value
    # https://dash.plot.ly/sharing-data-between-callbacks
    html.Div(id='json_structure', style={'display': 'none'})
])

modal_help = html.Div([
        dbc.Button("Help", className="btn-sm", id="open_modal"),
        dbc.Modal(
            [dbc.ModalHeader("How to upload a Structure"),
             dbc.ModalBody(dcc.Markdown("""
Use this interface to upload a file containing structural information.
Then select one of the other `actions` in the menu bar to operate on this structure.

The file extension (e.g. `.cif`) defines the format used to interpret the data.
All the formats supported by AbiPy and pymatgen are supported.
This includes Abinit netcdf files with the `.nc` extension, DDB files with the `_DDB` extension,
as well as Abinit input files with the `.abi` extension and Abinit output files with the `.abo` extension.
""")),
             dbc.ModalFooter(
                dbc.Button("Close", id="close", className="ml-auto")
             ),
            ], id="modal_help",
        ),
])


@app.callback(
    Output("modal_help", "is_open"),
    [Input("open_modal", "n_clicks"), Input("close", "n_clicks")],
    [State("modal_help", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


structure_upload = html.Div([
    dcc.Upload(id='upload-data',
        children=html.Div(['Drag and Drop or ', html.A('Select Files')]),
        style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Don't allow multiple files to be uploaded
        multiple=False,
    ),
    modal_help,
])


init_structure_layout = dbc.Container([
    dbc.Row([dbc.Col([html.H3("Upload file:"), structure_upload], md=4),
             dbc.Col([html.Div(id='print_structure')]),
            ]),
    ], className="mt-4",
)


input_groups = html.Div([
    dbc.InputGroup([
        dbc.InputGroupAddon("SpinMode:", addon_type="prepend"),
        dbc.Select(id="spin_mode", options=[
                {'label': "unpolarized", 'value': 'unpolarized'},
                {'label': "polarized", 'value': 'polarized'},
                {"label": "anti-ferromagnetic", "value": "afm"},
                {"label": "non-collinear with magnetism", "value": "spinor"},
                {"label": "non-collinear, no magnetism", "value": "spinor_nomag"},
            ], value="unpolarized"),
        ]),
    dbc.InputGroup([
        dbc.InputGroupAddon("Kppra:", addon_type="prepend"),
        dcc.Input(id="kppra", type="number", min=0, step=1000, value=1000),
        ]),
    dcc.RadioItems(id="gs_type",
        options=[
            {'label': 'SCF', 'value': 'scf'},
            {'label': 'Relax', 'value': 'relax'},
        ], value='scf'
    )
])

#'nosmearing': 1,
#'fermi_dirac': 3,
#'marzari4': 4,
#'marzari5': 5,
#'methfessel': 6,
#'gaussian': 7}

input_page_layout = html.Div([
    dbc.Row([dbc.Col(input_groups), dbc.Col(html.Div(id='gs_input_page_output'))]),
])

convert_format_input = dbc.InputGroup([
    dbc.InputGroupAddon("Output Format:", addon_type="prepend"),
    dbc.Select(id="convert_format", options=[
                    {"label": k, "value": k} for k in "abivars,cif,xsf,poscar,qe,siesta,wannier90,cssr,json".split(",")
                    ], value="cif"),
])

convert_page_layout = html.Div([
    dbc.Row([dbc.Col(convert_format_input),
    dbc.Col(html.Div(id='convert_page_output'))]),
])

symmetry_options_input = dbc.InputGroup([
    dbc.InputGroupAddon("spglib symprec (A):", addon_type="prepend"),
    dbc.Input(id="spglib_symprec", type="number", min=1e-6, step=0.02, value=0.01),
    dbc.InputGroupAddon("spglib angdeg (degrees):", addon_type="prepend"),
    dbc.Input(id="spglib_angtol", type="number", min=1e-6, step=1, value=5),
    #dbc.InputGroupAddon("Abinit tolsym:", addon_type="prepend"),
    #dbc.Input(id="abinit_tolsym", type="number", min=1e-12, step=1e-1, value=1e-8),
])


symmetry_analysis_page_layout = html.Div([
    dbc.Row(dbc.Col(symmetry_options_input)),
    dbc.Row(html.Div(id='symmetry_page_output')),
])


@app.callback(Output('gs_input_page_output', 'children'),
              [Input('spin_mode', 'value'),
               Input('kppra', 'value'),
               Input('gs_type', 'value'),
              ],
              [State("json_structure", "children")],
              )
def update_gs_input(spin_mode, kppra, gs_type, json_structure):
    if not json_structure: return structure_undefined_error()

    structure = abilab.mjson_loads(json_structure)
    pseudos = os.path.join(abidata.pseudo_dir, "14si.pspnc")

    # Build input for GS calculation
    # ecut must be specified because this pseudopotential does not provide hints for ecut.
    try:
        gs_inp = abilab.gs_input(
            structure, pseudos,
            kppa=kppra, ecut=8, spin_mode=spin_mode, smearing=None)
        gs_inp.pop_vars(("charge", "chksymbreak"))

        if gs_type == "relax":
            gs_inp.set_vars(optcell=2, ionmov=2, ecutsm=0.5, dilatmx=1.05)

        #multi = ebands_input(structure, pseudos,
        #                 kppa=kppra, nscf_nband=None, ndivsm=15,
        #                 ecut=8, pawecutdg=None, scf_nband=None, accuracy="normal", spin_mode=spin_mode,
        #                 smearing="fermi_dirac:0.1 eV", charge=None, dos_kppa=None):

        gs_inp.set_mnemonics(False)
    except Exception as exc:
        return html.Div([dbc.Jumbotron([html.H2("There was an error processing this file. %s" % str(exc),
                         className="text-danger")])])

    s = DangerouslySetInnerHTML(gs_inp._repr_html_())

    return html.Div([
        copy_to_clipboard_button(id="copy_gsinput"),
        html.Hr(),
        html.Pre(s, id="copy_gsinput"), # className="text-sm-left",
        ])


@app.callback(Output('convert_page_output', 'children'),
              [Input('convert_format', 'value')],
              [State("json_structure", "children")],
             )
def update_convert_page(output_format, json_structure):
    if not json_structure: return structure_undefined_error()
    structure = abilab.mjson_loads(json_structure)
    return html.Div([
        copy_to_clipboard_button(id="copy_convert_output"),
        html.Hr(),
        html.Pre(structure.convert(fmt=output_format), id="copy_convert_output"),
    ])


@app.callback(Output('symmetry_page_output', 'children'),
              [Input('spglib_symprec', 'value'),
               Input('spglib_angtol', 'value'),
               #Input('abinit_tolsym', 'value'),
              ],
              [State("json_structure", "children")],
             )
def update_symmetry_page(spglib_symprec, spglib_angtol, json_structure):
    if not json_structure: return structure_undefined_error()
    structure = abilab.mjson_loads(json_structure)

    #print(f"symprec: {spglib_symprec}, spglib_angtol: {spglib_angtol}, abinit_tolsym: {abinit_tolsym}")
    s = structure.spget_summary(symprec=spglib_symprec, angle_tolerance=spglib_angtol)
    return html.Div([html.Br(), html.Pre(s)])


kpath_options_input = dbc.InputGroup([
    dbc.InputGroupAddon("Output Format:", addon_type="prepend"),
    dbc.Select(id="kpath_fmt",
               options=[{"label": k, "value": k} for k in ("abinit", "wannier90", "siesta")], value="abinit"),
    dbc.InputGroupAddon("Line Density:", addon_type="prepend"),
    dbc.Input(id="kpath_line_density", min=1, step=10, value=10, type="number"),
])


kpath_page_layout = html.Div([
    dbc.Row(dbc.Col(html.Div(kpath_options_input))),
    dbc.Row(html.Div(id='kpath_page_output')),
])


@app.callback(Output('kpath_page_output', 'children'),
              [Input('kpath_fmt', 'value'),
               Input('kpath_line_density', 'value'),
              ],
              [State("json_structure", "children")],
             )
def update_kpath_page(kpath_fmt, kpath_line_density, json_structure):
    if not json_structure: return structure_undefined_error()
    structure = abilab.mjson_loads(json_structure)

    #print(f"kpath_fmt: {kpath_fmt}, kpath_line_density: {kpath_line_density}")
    try:
        s = structure.get_kpath_input_string(fmt=kpath_fmt, line_density=int(kpath_line_density))
    except Exception as exc:
        return div_with_exc(exc)

    return html.Div([
        html.Br(),
        copy_to_clipboard_button(id="kpath_string"),
        html.Hr(),
        html.Pre(s, id="kpath_string")
        ])


def structure_undefined_error():
    return html.Div([
        html.Br(),
        dbc.Jumbotron([html.H2("Use `Select Structure` to upload the file", className="text-danger")])
    ])


@app.callback([Output('print_structure', 'children'),
               Output('json_structure', 'children')],
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified')])
def structure_from_upload(contents, filename, date):
    #print("filename, date", filename, date)
    content_type, content_string = contents.split(',')

    if filename.endswith(".nc"):
        mode = "wb"
        decoded = base64.b64decode(content_string)
    else:
        mode = "wt"
        decoded = base64.b64decode(content_string).decode('utf-8')

    try:
        with tempfile.TemporaryDirectory() as tmp_dirname:
            tmp_path = os.path.join(tmp_dirname, filename)
            with open(tmp_path, mode) as fh:
                fh.write(decoded)
            structure = abilab.Structure.from_file(tmp_path)

        return html.Div([html.Br(), html.Pre(str(structure))]), structure.to_json()

    except Exception as exc:
        return html.Div([dbc.Jumbotron([html.H2("There was an error processing this file. %s" % str(exc),
                         className="text-danger")])]), None


about_page_layout = html.Div([
    dcc.Markdown("""
This web app is built on the
[AbiPy](https://github.com/abinit/abipy),
[pymatgen](https://github.com/materialsproject/pymatgen) and
[spglib](https://atztogo.github.io/spglib/) libraries.

Note that this tool provides a simplified interface to the AbiPy API.
For a more flexible interface, please use the AbiPy objects to generate input files and workflows.

**Privacy**: We do not store any files or data on our server.
Any uploaded files is analyzed and the results presented immediately.
""")
])

@app.callback(Output("page-content", "children"),
              [Input("url", "pathname")])
def render_page_content(pathname):

    if pathname in ["/", "/page-1"]:
        return init_structure_layout
    elif pathname == "/page-2":
        return convert_page_layout
    elif pathname == "/page-3":
        return symmetry_analysis_page_layout
    elif pathname == "/page-4":
        return kpath_page_layout
    elif pathname == "/page-5":
        return input_page_layout
    elif pathname == "/page-6":
        return about_page_layout

    # If the user tries to reach a different page, return a 404 message
    return dbc.Jumbotron([html.H1("404: Not found", className="text-danger"),
                          html.Hr(),
                          html.P(f"The pathname {pathname} was not recognised..."),
                         ])


# this callback uses the current pathname to set the active state of the
# corresponding nav link to true, allowing users to tell see page they are on
@app.callback(
    [Output(f"page-{i}-link", "active") for i in range(1, num_pages + 1)],
    [Input("url", "pathname")],
)
def toggle_active_links(pathname):
    if pathname == "/":
        # Treat page 1 as the homepage / index
        out = num_pages * [False]
        out[0] = True
        return out
    return [pathname == f"/page-{i}" for i in range(1, num_pages + 1)]


@app.callback(
    Output("sidebar", "className"),
    [Input("sidebar-toggle", "n_clicks")],
    [State("sidebar", "className")],
)
def toggle_classname(n, classname):
    if n and classname == "":
        return "collapsed"
    return ""


@app.callback(
    Output("collapse", "is_open"),
    [Input("navbar-toggle", "n_clicks")],
    [State("collapse", "is_open")],
)
def toggle_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


if __name__ == "__main__":
    app.run_server(debug=True)
