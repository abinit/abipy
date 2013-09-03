from __future__ import division, print_function

from pymatgen.io.abinitio.workflow import Workflow as pmWorkflow

class Workflow(pmWorkflow):

    def wxshow_inputs(self):
        """Open a noteboox dysplaying the input files of the workflow."""
        from abipy.gui.wxapps import wxapp_showfiles
        wxapp_showfiles(dirpath=self.workdir, walk=True, wildcard="*.abi").MainLoop()

    def wxshow_outputs(self):
        """Open a noteboox dysplaying the output files of the workflow."""
        from abipy.gui.wxapps import wxapp_showfiles
        wxapp_showfiles(dirpath=work.workdir, walk=True, wildcard="*.abo").MainLoop()

    #def wxshow_events(self):
    #    """Open a noteboox dysplaying the events (warnings, comments, warnings) of the workflow."""
    #    from abipy.gui.wxapps import wxapp_showevents
    #    wxapp_showfiles(dirpath=work.workdir, walk=True, wildcard="*.abo").MainLoop()
