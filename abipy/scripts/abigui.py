#!/usr/bin/env python
"""
Script to start the panel-based AbiPy web GUI.
"""
import os
import click
import panel as pn


@click.command()
@click.option('--port', default=5006, help="Port to listen on. Default: 5006")
@click.option('--address', default=None, help="The address the server should listen on for HTTP requests.")
@click.option('--show', default=True, help="Open server app(s) in a browser")
@click.option('--num_procs', default=1, help="Number of worker processes for the app. Defaults to 1")
@click.option('--panel-template', default="FastList",
              help="Specify template for panel dasboard." +
                   "Possible values are: FastList, FastGrid, Golden, Bootstrap, Material, React, Vanilla." +
                   "Default: FastList")
def gui_app(port, address, show, num_procs, panel_template):

    from abipy.panels.core import abipanel, get_abinit_template_cls_kwds #, AbipyParameterized
    import abipy.panels as mod

    # Load abipy/panel extensions.
    abipanel(panel_template=panel_template)

    assets_path = os.path.join(os.path.dirname(mod.__file__), "assets")
    #AbipyParameterized.uses_abinit_server = True # TODO

    main_home = pn.Column(pn.pane.Markdown("""

![AbiPy Logo](assets/img/abipy_logo.png)

# AbiPy Web App

This web application exposes some of the capabilities of the AbiPy package.
It consists of multiple pages each of which provides specialized tools to operate on a particular ABINIT file.
To access one of these tools, click one of the links in the sidebar or, alternatively, use the links below.
To open/close the sidebar, click on the Hamburger Menu Icon ☰ in the header.

## [Abinit Input Generator](/input_generator)

Generate ABINIT input files for performing basic

   - ground-state calculations
   - band structure calculations
   - DFPT phonon calculations

starting from a crystalline structure provided by the user either via an external file
or through the Materials Project identifier (*mp-id*)

## [Structure Analyzer](/structure_analyzer)

This application allows user to upload a file with structural info and operate on it.

## [Output File Analyzer](/outfile)

Post-process the data stored in one of the ABINIT output files.

## [DDB File Analyzer](/ddb)

This application allows users to post-process the data stored in one of the Abinit output files.
The main difference with respect to [Abinit Output File Analyzer](/outfile) is that
it is also possible to fetch the DDB file from the Materials Project Database.

## [Abo File Analyzer](/abo)

Analyze the Abinit main output file

## [Ebands vs MP](/ebands_vs_mp)

Compare your electronic bands with the MP

## [DDB vs MP](/ddb_vs_mp)

Compare your DDB with the MP

## [SKW Analyzer](/skw)

This tool allows one to interpolate the KS energies with the star-function method and
compare the interpolated band structure with the *ab-initio* one.

""", sizing_mode="stretch_width"),
    sizing_mode="stretch_width")

    cls, kwds = get_abinit_template_cls_kwds()
    print("Using panel template:", cls)
    home = cls(main=main_home, title="AbiPy GUI Home", **kwds)

    # Import the apps and define routies for each page.
    from abipy.panels.structure import InputFileGenerator
    from abipy.panels.ddb import PanelWithFileInput, DdbPanelWithFileInput, CompareDdbWithMP
    from abipy.panels.electrons import SkwPanelWithFileInput, CompareEbandsWithMP
    from abipy.panels.outputs import AbinitOutputFilePanelWithFileInput as abo_cls

    # url --> (app, title)
    app_routes_titles = {
        "/": (home, "Home"),
        "/input_generator": (InputFileGenerator().get_panel(), "Abinit Input Generator"),
        "/structure_analyzer": (PanelWithFileInput(use_structure=True).get_panel(), "Structure Analyzer"),
        "/outfile": (PanelWithFileInput().get_panel(), "Output File Analyzer"),
        "/ddb": (DdbPanelWithFileInput().get_panel(), "DDB File Analyzer"),
        "/abo": (abo_cls().get_panel(), "Abo File Analyzer"),
        "/ebands_vs_mp": (CompareEbandsWithMP().get_panel(), "Compare Ebands with MP"),
        "/ddb_vs_mp": (CompareDdbWithMP().get_panel(), "Compare DDB with MP"),
        "/skw": (SkwPanelWithFileInput().get_panel(), "SKW Analyzer"),
        #"/abilog": (PanelWithFileInput().get_panel(), "DDB File Analyzer"),
        #"/gs_autoparal": (PanelWithFileInput().get_panel(), "DDB File Analyzer"),
    }

    app_routes = {k: v[0] for (k, v) in app_routes_titles.items()}
    app_title = {k: v[1] for (k, v) in app_routes_titles.items()}

    # Add links to sidebar of each app so that we can navigate easily.
    links = "\n".join(f"- [{title}]({url})" for (url, title) in app_title.items())
    links = pn.pane.Markdown(links)

    for url, app in app_routes.items():
        #if not hasattr(app, "sidebar"): continue
        app.sidebar.append(links)

    # Call pn.serve to serve the multipage app.
    kwargs = dict(address=address,
                  port=port,
                  #dev=True,
                  #start=True,
                  show=show,
                  #title=app_title,
                  num_procs=num_procs,
                  static_dirs={"/assets": assets_path},
    )
    from pprint import pformat
    print("Calling pn.serve with kwargs:\n", pformat(kwargs))

    pn.serve(app_routes, **kwargs)


if __name__ == "__main__":
    gui_app()
