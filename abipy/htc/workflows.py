from __future__ import division, print_function

from pymatgen.io.abinitio.workflow import Workflow as pmWorkflow

class Workflow(pmWorkflow):

    def wxshow_inputs(self):
        from abipy.gui.wxapps import wxapp_showfiles
        wxapp_showfiles(dirpath=self.workdir, walk=True, wildcard="*.abi").MainLoop()

    def wxshow_outputs(self):
        from abipy.gui.wxapps import wxapp_showfiles
        wxapp_showfiles(dirpath=work.workdir, walk=True, wildcard="*.abo").MainLoop()
