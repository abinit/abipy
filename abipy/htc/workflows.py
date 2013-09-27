from __future__ import division, print_function

from pymatgen.io.abinitio.workflow import IterativeWorkflow, Workflow as pmWorkflow
from abipy.htc.abitimer import AbinitTimerParser

class Workflow(pmWorkflow):
    """Hook used to add and test additional features to the pymatgen workflow."""

    def parse_timers(self):
        """Parse the TIMER section reported in the ABINIT output files."""
        filenames = [task.output_file.path for task in self]

        parser = AbinitTimerParser()
        parser.parse(filenames)

        return parser
