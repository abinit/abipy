import os
import platform
import panel as pn
from abipy import abilab

abilab.abipanel()
#from abipy.panel.core import AbipyParameterized
#AbipyParameterized.uses_abinit_server = True

#logos = pn.Row(
#    #pn.pane.PNG("assets/abinit_logo.png", width=500, link_url="https://github.com/abinit/abinit"),
#    pn.pane.PNG("assets/abipy_logo.png", width=500, link_url="https://github.com/abinit/abipy"),
#    #ncols=2,
#)

alert = pn.pane.Alert("""
Please note that this web interface is not designed to handle large data transfer.
To post-process the data stored in a big file e.g. a WFK.nc file,
we strongly suggest running AbiPy on the same machine where the file is hosted.

Examples of post-processing scripts are available in the
[AbiPy gallery](https://abinit.github.io/abipy/gallery/index.html).
""", alert_type="danger")


main_home = pn.Column(
    #logos,
    #pn.pane.PNG("assets/abipy_logo.png", width=500, link_url="https://github.com/abinit/abipy"),
    pn.pane.Markdown("""

![Caption](assets/abipy_logo.png)

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

## [SKW interpolation](/skw)

This tool allows one to interpolate the KS energies with the star-function method and
compare the interpolated band structure with the *ab-initio* one.

""", sizing_mode="stretch_width"),
alert,
sizing_mode="stretch_width")


home = pn.template.FastListTemplate(main=main_home, title="AbiPy GUI Home",
                                    header_background="#ff8c00") # Dark orange

from abipy.panels.structure import InputFileGenerator
from abipy.panels.ddb import PanelWithFileInput, DdbPanelWithFileInput, CompareDdbWithMP
from abipy.panels.electrons import SkwPanelWithFileInput, CompareEbandsWithMP
from abipy.panels.outputs import AbinitOutputFilePanelWithFileInput as abo_cls

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

links = "\n".join(f"- [{title}]({url})" for (url, title) in app_title.items())
#print(links)
links = pn.pane.Markdown(links)

for url, app in app_routes.items():
    app.sidebar.append(links)

kwargs = dict(#address=os.getenv("BOKEH_ADDRESS", "0.0.0.0"),
              port=80,
              dev=True,
              start=True,
              show=True,
              title=app_title,
              #num_procs=4 if platform.system() != "Windows" else 1,
              static_dirs={"/assets": "assets"},
              )

pn.serve(app_routes, **kwargs)
