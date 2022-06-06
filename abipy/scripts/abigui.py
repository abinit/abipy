#!/usr/bin/env python
"""
Script to start the panel-based AbiPy web GUI.
"""
import sys
import os
import argparse
import panel as pn
import abipy.tools.cli_parsers as cli

from pprint import pformat
from abipy.core.release import version

def main():

    def str_examples():
        return """\
Usage example:

    abigui.py --verbose
    abigui.py --num_procs 4 --has-remote-server  => Options used for gui.abinit.org

NB: To deploy the GUI, use the ~abipy/dev_scripts/deploy_abigui.sh script.
"""

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    def get_copts_parser():

        # Parent parser implementing common options.
        p = argparse.ArgumentParser(add_help=False)
        p.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                       help='Verbose, can be supplied multiple times to increase verbosity')

        p.add_argument('--loglevel', default="ERROR", type=str,
                            help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

        from abipy.core.release import __version__
        p.add_argument('-V', '--version', action='version', version=__version__)

        return p

    copts_parser = get_copts_parser()

    #cli.customize_mpl(options)

    parents = [copts_parser, cli.pn_serve_parser()] #, plot_parser]

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=str_examples(), parents=parents,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    from abipy.panels.core import abipanel, get_abinit_template_cls_kwds, AbipyParameterized

    # Load abipy/panel extensions and set the default template
    #tmpl_kwds.update(dict(
    #    sidebar_width=240,
    #    #sidebar_width=280,
    #    #background_color="yellow",
    #))
    abipanel(panel_template=options.panel_template)

    if options.has_remote_server:
        print("Enforcing limitations on what the user can do on the abinit server")
        AbipyParameterized.has_remote_server = options.has_remote_server

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

To access the tools, click one of the links in the sidebar or, alternatively, use the links below.

To **open/close** the sidebar, click on the Hamburger Menu Icon ☰ in the header.

Note that the **file extension matters** as the GUI won't work properly if you upload files
with extensions that are not recognized by AbiPy.

"""
    tmpl_cls, tmpl_kwds = get_abinit_template_cls_kwds()

    if options.verbose:
        print("Using default template:", tmpl_cls, "with kwds:\n", pformat(tmpl_kwds), "\n")

    # url --> (cls, title)
    app_routes_titles = {
        "/": (tmpl_cls, "AbiPy GUI Home"),
        "/input_generator": (InputFileGenerator, "Abinit Input Generator"),
        "/structure_analyzer": (PanelWithStructureInput, "Structure Analyzer"),
        "/outfile": (PanelWithFileInput, "Output File Analyzer"),
        "/ddb": (DdbPanelWithFileInput, "DDB File Analyzer"),
        "/abo": (abo_cls, "Abo File Analyzer"),
        "/ebands_vs_mp": (CompareEbandsWithMP, "Compare Ebands with MP"),
        "/ddb_vs_mp": (CompareDdbWithMP, "Compare DDB with MP"),
        #"/abilog": (PanelWithFileInput().get_panel(), "DDB File Analyzer"),
        #"/state": (pn.state, "State"),
    }

    if not options.has_remote_server:
        # Add additional apps.
        app_routes_titles.update({
            "/skw": (SkwPanelWithFileInput, "SKW Analyzer"),
            "/robot": (RobotWithFileInput, "Robot Analyzer"),
        })

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

    class AppBuilder:

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
    serve_kwargs = cli.get_pn_serve_kwargs(options)

    pn.serve(app_routes, **serve_kwargs)


if __name__ == "__main__":
    sys.exit(main())
