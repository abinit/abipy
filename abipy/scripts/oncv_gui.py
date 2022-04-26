#!/usr/bin/env python
"""
Script to start the panel-based Oncvpsp web GUI.
"""
import sys
import os
import click
import panel as pn

from pprint import pformat
from abipy.core.release import version
from abipy.panels.oncvpsp_gui import OncvGui


@click.command()
@click.argument('filepath')
@click.option('--port', default=5006, help="Port to listen on. Default: 5006")
@click.option('--address', default=None, help="The address the server should listen on for HTTP requests.")
@click.option('--show', default=True, help="Open server app(s) in the web browser")
@click.option('--num_procs', default=1, help="Number of worker processes for the app. Defaults to 1")
@click.option('--panel-template', "-pnt", default="FastList",
              help="Specify template for panel dasboard." +
                   "Possible values are: FastList, FastGrid, Golden, Bootstrap, Material, React, Vanilla." +
                   "Default: FastList")
@click.option('--has-remote-server', default=False, is_flag=True,
              help="True if we are running on the ABINIT server. This flag activates limitations on what the user can do." +
                   "Default: False")
@click.option("--websocket-origin", default=None, type=str,
        help="Public hostnames which may connect to the Bokeh websocket.\n Syntax: " +
              "HOST[:PORT] or *. Default: None")
@click.option('--max_size_mb', default=150, type=int,
        help="Maximum message size in Mb allowed by Bokeh and Tornado. Default: 150")
@click.option("-v", '--verbose', default=0, count=True, help="Verbosity level")
@click.version_option(version=version, message='%(version)s')
def gui_app(filepath, port, address, show, num_procs, panel_template, has_remote_server, websocket_origin, max_size_mb, verbose):

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

    tmpl_cls, tmpl_kwds = get_abinit_template_cls_kwds()

    if verbose:
        print("Using default template:", tmpl_cls, "with kwds:\n", pformat(tmpl_kwds), "\n")

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

    #if options.seaborn:
    # Use seaborn settings.
    import seaborn as sns
    context = "paper"
    context = "notebook"
    sns.set(context=context, style='darkgrid', palette='deep',
            font='sans-serif', font_scale=1, color_codes=False, rc=None)

    #pn.extension('ace')
    #pn.extension('ipywidgets')

    def build():
        gui = OncvGui.from_file(os.path.abspath(filepath))
        template = tmpl_cls(main=gui.get_panel(), title="Oncvpsp GUI", **tmpl_kwds)
        return template

    #if verbose:
    print("Calling pn.serve with serve_kwargs:\n", pformat(serve_kwargs), "\n")

    return pn.serve(build, **serve_kwargs)


if __name__ == "__main__":
    gui_app()
