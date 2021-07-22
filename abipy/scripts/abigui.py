#!/usr/bin/env python
"""
Script to start the panel-based AbiPy web GUI.
"""
import os
import click
import panel as pn

from abipy.core.release import version


@click.command()
@click.option('--port', default=5006, help="Port to listen on. Default: 5006")
@click.option('--address', default=None, help="The address the server should listen on for HTTP requests.")
@click.option('--show', default=True, help="Open server app(s) in a browser")
@click.option('--num_procs', default=1, help="Number of worker processes for the app. Defaults to 1")
@click.option('--panel-template', "-pnt", default="FastList",
              help="Specify template for panel dasboard." +
                   "Possible values are: FastList, FastGrid, Golden, Bootstrap, Material, React, Vanilla." +
                   "Default: FastList")
@click.option('--has-remote-server', default=False, is_flag=True,
              help="True if we are running on the ABINIT server. This flag activates limitations on what the user can do." +
                   "Default: False")
@click.version_option(version=version, message='%(version)s')
def gui_app(port, address, show, num_procs, panel_template, has_remote_server):

    from abipy.panels.core import abipanel, get_abinit_template_cls_kwds, AbipyParameterized
    import abipy.panels as mod

    # Load abipy/panel extensions.
    abipanel(panel_template=panel_template)

    assets_path = os.path.join(os.path.dirname(mod.__file__), "assets")
    print("has_remote_server:", has_remote_server)
    if has_remote_server:
        print("Enforce limitations on what the user can do on the abinit server")
        AbipyParameterized.has_remote_server = has_remote_server # TODO

    # Import the apps and define routies for each page.
    from abipy.panels.structure import InputFileGenerator
    from abipy.panels.ddb import (PanelWithFileInput, PanelWithStructureInput, DdbPanelWithFileInput, CompareDdbWithMP,
                                  RobotWithFileInput)
    from abipy.panels.electrons import SkwPanelWithFileInput, CompareEbandsWithMP
    from abipy.panels.outputs import AbinitOutputFilePanelWithFileInput as abo_cls

    intro = """
![AbiPy Logo](assets/img/abipy_logo.png)

# AbiPy Web App

This web application exposes some of the capabilities of the [AbiPy package](https://github.com/abinit/abipy)
It consists of **multiple pages** each of which provides **specialized tools** to operate on a particular ABINIT file.
To access one of these tools, click one of the links in the sidebar or, alternatively, use the links below.
To open/close the sidebar, click on the Hamburger Menu Icon ☰ in the header.

Note that the **file extension** matters as the GUI won't work properly if you upload files
with extensions that are not recognized by AbiPy.

"""
    cls, cls_kwds = get_abinit_template_cls_kwds()
    print("Using panel template:", cls)

    # url --> (cls, title)
    app_routes_titles = {
        "/": (cls, "AbiPy GUI Home"),
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
        #"/gs_autoparal": (PanelWithFileInput().get_panel(), "DDB File Analyzer"),
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

    main_home = pn.Column(pn.pane.Markdown(intro, width=600, sizing_mode="stretch_width"),
                          sizing_mode="stretch_both")

    # Add links to sidebar of each app so that we can navigate easily.
    links = "\n".join(f"- [{title}]({url})" for (url, title) in app_title.items())
    links = pn.Column(pn.pane.Markdown(links))

    def func(cls, **cls_kwargs):
        app = cls(**cls_kwargs)
        if hasattr(app, "get_panel"):
            app = app.get_panel()
        app.sidebar.append(links)
        #app.header.append(links)
        return app

    class Partial():
        # https://stackoverflow.com/questions/45485017/why-does-functools-partial-not-detected-as-a-types-functiontype
        def __init__(self, func, *args, **kwargs):
            self.func = func
            self.args = args
            self.kwargs = kwargs

        def view(self):
            return self.func(*self.args, **self.kwargs)


    cls_kwds.update(dict(
        sidebar_width=240,
        #sidebar_width=280,
        #background_color="yellow",
    ))


    for url, cls in app_routes.items():
        if url == "/":
            app_routes[url] = Partial(func, cls, main=main_home, **cls_kwds).view
        elif url == "/state":
            app_routes[url] = cls
        else:
            app_routes[url] = Partial(func, cls, **cls_kwds).view

    # Call pn.serve to serve the multipage app.
    serve_kwargs = dict(address=address,
                  port=port,
                  #dev=True,
                  #start=True,
                  show=show,
                  #title=app_title,
                  num_procs=num_procs,
                  static_dirs={"/assets": assets_path},
    )

    from pprint import pformat
    print("Calling pn.serve with serve_kwargs:\n", pformat(serve_kwargs))

    pn.serve(app_routes, **serve_kwargs)


if __name__ == "__main__":
    gui_app()
