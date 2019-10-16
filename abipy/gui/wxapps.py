"""This module import all WxWidgets applications provides by Abipy in a single namespace."""

from .events import wxapp_events
from .browser import wxapp_dirbrowser, wxapp_listbrowser
from .scissors import wxapp_scissors
from .wfkviewer import wxapp_wfkviewer
from .sigresviewer import wxapp_sigresviewer
from .comparison import wxapp_comparison
from .editor import wxapp_showfiles
from .converter import wxapp_converter
from .structure import wxapp_structure_converter
from .fftprof import wxapp_fftprof
from .flowviewer import wxapp_flow_viewer
from .gsrviewer import wxapp_gsrviewer
from .mdfviewer import wxapp_mdfviewer
from .wxncview import wxapp_ncview
#from .oncvgui import wxapp_oncvpsp
from .phbstviewer import wxapp_phbstviewer


# Map abinit file extensions to WX Applications.
_EXT2APP = {
    "WFK-etsf.nc": wxapp_wfkviewer,
    "SIGRES.nc": wxapp_sigresviewer,
    "GSR.nc": wxapp_gsrviewer,
    "MDF.nc": wxapp_mdfviewer,
    "PHBST.nc": wxapp_phbstviewer,
    #".abi": MyEditorApp,
    #".abo": MyEditorApp,
    #".log": MyEditorApp,
    #".sh": MyEditorApp,
    #".err": MyEditorApp,
    #".files": MyEditorApp,
}


def file2appcls(filepath):
    import os
    ext = filepath.split("_")[-1]
    try:
        return _EXT2APP[ext]

    except KeyError:
        root, ext = os.path.splitext(filepath)
        try:
            return _EXT2APP[ext]
        except KeyError:
            # No frame registered for the file.
            # Open NcViewer if we have a netcdf file else None
            if filepath.endswith(".nc"): return wxapp_ncview
            return None
