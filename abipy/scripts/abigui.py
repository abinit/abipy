#!/usr/bin/env python
"""
Script to start the panel-based AbiPy web GUI.
"""
import os
import click
import panel as pn

from pprint import pformat
from abipy.core.release import version


@click.command()
@click.option('--port', default=5006, help="Port to listen on. Default: 5006")
@click.option('--address', default=None, help="The address the server should listen on for HTTP requests.")
@click.option('--show', default=True, help="Open server app(s) in the web browser")
@click.option('--num_procs', default=1, help="Number of worker processes for the app. Defaults to 1")
@click.option('--panel-template', "-pnt", default="FastList",
              help="Specify template for panel dasboard." +
                   "Possible values are: FastList, FastGrid, Golden, Bootstrap, Material, React, Vanilla. " +
                   "Default: FastList")
@click.option('--has-remote-server', default=False, is_flag=True,
              help="True if we are running on the ABINIT server. This flag activates limitations on what the user can do. " +
                   "Default: False")
@click.option("--websocket-origin", default=None, type=str,
        help="Public hostnames which may connect to the Bokeh websocket.\n Syntax: " +
              "HOST[:PORT] or *. Default: None")
@click.option('--max_size_mb', default=150, type=int,
        help="Maximum message size in Mb allowed by Bokeh and Tornado. Default: 150")
@click.option("-v", '--verbose', default=0, count=True, help="Verbosity level")
@click.version_option(version=version, message='%(version)s')
def gui_app(port, address, show, num_procs, panel_template, has_remote_server, websocket_origin, max_size_mb, verbose):

    from abipy.panels.core import abipanel, get_abinit_template_cls_kwds, AbipyParameterized
    import abipy.panels as mod
    assets_path = os.path.join(os.path.dirname(mod.__file__), "assets")

    # Load abipy/panel extensions and set the default template
    #tmpl_kwds.update(dict(
    #    sidebar_width=240,
    #    #sidebar_width=280,
    #    #background_color="yellow",
    #))
    abipanel(panel_template=panel_template)

    if has_remote_server:
        print("has_remote_server:", has_remote_server)
        print("Enforcing limitations on what the user can do on the abinit server")
        # TODO finalize, remove files created by user
        AbipyParameterized.has_remote_server = has_remote_server

    # Import the apps and define routes for each page.
    from abipy.panels.structure import InputFileGenerator
    from abipy.panels.ddb import (PanelWithFileInput, PanelWithStructureInput, DdbPanelWithFileInput, CompareDdbWithMP,
                                  RobotWithFileInput)
    from abipy.panels.electrons import SkwPanelWithFileInput, CompareEbandsWithMP
    from abipy.panels.outputs import AbinitOutputFilePanelWithFileInput as abo_cls

    intro = """
![AbiPy Logo](assets/img/abipy_logo.png)

# AbiPy Web App

This web application exposes some of the capabilities of the [AbiPy package](https://github.com/abinit/abipy).
It consists of **multiple pages** each of which provides **specialized tools** to operate on a particular ABINIT file.
To access one of these tools, click one of the links in the sidebar or, alternatively, use the links below.
To **open/close** the sidebar, click on the Hamburger Menu Icon ☰ in the header.

Note that the **file extension** matters as the GUI won't work properly if you upload files
with extensions that are not recognized by AbiPy.

"""
    tmpl_cls, tmpl_kwds = get_abinit_template_cls_kwds()

    if verbose:
        print("Using default template:", tmpl_cls, "with kwds:\n", pformat(tmpl_kwds), "\n")

    # url --> (cls, title)
    app_routes_titles = {
        "/": (tmpl_cls, "AbiPy GUI Home"),
        "/input_generator": (InputFileGenerator, "Abinit Input Generator"),
        #"/structure_analyzer": (PanelWithFileInput(use_structure=True).get_panel(), "Structure Analyzer"),
        "/structure_analyzer": (PanelWithStructureInput, "Structure Analyzer"),
        "/outfile": (PanelWithFileInput, "Output File Analyzer"),
        "/ddb": (DdbPanelWithFileInput, "DDB File Analyzer"),
        "/abo": (abo_cls, "Abo File Analyzer"),
        "/ebands_vs_mp": (CompareEbandsWithMP, "Compare Ebands with MP"),
        "/ddb_vs_mp": (CompareDdbWithMP, "Compare DDB with MP"),
        "/skw": (SkwPanelWithFileInput, "SKW Analyzer"),
        "/robot": (RobotWithFileInput, "Robot Analyzer"),
        #"/abilog": (PanelWithFileInput().get_panel(), "DDB File Analyzer"),
        #"/state": (pn.state, "State"),
    }

    app_routes = {k: v[0] for (k, v) in app_routes_titles.items()}
    app_title = {k: v[1] for (k, v) in app_routes_titles.items()}

    for url, (cls, title) in app_routes_titles.items():
        if url in ("/", "/state"): continue
        intro += f"""

### [{title}]({url})

{cls.info_str}
"""

    main_home = pn.Column(pn.pane.Markdown(intro, sizing_mode="stretch_both"),
                          sizing_mode="stretch_both")

    # Add links to sidebar of each app so that we can navigate easily.
    links = "\n".join(f"- [{title}]({url})" for (url, title) in app_title.items())
    links = pn.Column(pn.pane.Markdown(links))

    class AppBuilder():

        def __init__(self, app_cls, sidebar_links, app_kwargs=None):
            self.app_cls = app_cls
            self.sidebar_links = sidebar_links
            self.app_kwargs = app_kwargs if app_kwargs else {}

        def build(self):
            app = self.app_cls(**self.app_kwargs)
            if hasattr(app, "get_panel"):
                app = app.get_panel()

            if hasattr(app, "sidebar") and self.sidebar_links:
                app.sidebar.append(self.sidebar_links)
                #app.header.append(self.sidebar_links)

            return app

    for url, app_cls in app_routes.items():
        if url == "/":
            app_kwargs = {"main": main_home}
            app_kwargs.update(tmpl_kwds)
            app_routes[url] = AppBuilder(app_cls, sidebar_links=links, app_kwargs=app_kwargs).build
        elif url == "/state":
            app_routes[url] = app_cls
        else:
            app_routes[url] = AppBuilder(app_cls, sidebar_links=links).build

    # Call pn.serve to serve the multipage app.
    serve_kwargs = dict(
        address=address,
        port=port,
        #dev=True,
        #start=True,
        show=show,
        debug=verbose > 0,
        #title=app_title,
        num_procs=num_procs,
        static_dirs={"/assets": assets_path},
        websocket_origin=websocket_origin,
        #websocket_origin="*",
        #
        # Increase the maximum websocket message size allowed by Bokeh
        # https://panel.holoviz.org/reference/widgets/FileInput.html
        websocket_max_message_size=max_size_mb * 1024**2,
        # Increase the maximum buffer size allowed by Tornado
        http_server_kwargs={'max_buffer_size': max_size_mb * 1024**2},
    )

    #if verbose:
    print("Calling pn.serve with serve_kwargs:\n", pformat(serve_kwargs), "\n")

    pn.serve(app_routes, **serve_kwargs)


if __name__ == "__main__":
    gui_app()
