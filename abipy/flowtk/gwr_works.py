# coding: utf-8
"""
Works and Flows for GWR calculations (GW with supercells).
"""
from __future__ import annotations

#import numpy as np

from abipy.abio.inputs import AbinitInput, RUNL, GWR_TASK
from abipy.electrons.gwr import GwrRobot
from .nodes import Node
from .tasks import TaskManager # ScfTask, NscfTask,
from .works import Work


class DirectDiagoWork(Work):
    """
    This work performs the direct diagonalization of the KS Hamiltonian
    using the density produced by a GS-SCF run and produces a WFK file with
    empty states in the outdir of the second task.

    NB: An Abinit build with Scalapack is required to run GWR.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: DirectDiagoWork
    """

    @classmethod
    def from_scf_input(cls, scf_input: AbinitInput, green_nband: int,
                       manager: TaskManager=None) -> DirectDiagoWork:
        """
        Args:
            scf_input: Input for the GS-SCF calculations.
            green_nband: Number of bands to compute in the direct diagonalization.
                A negative value activate full diagonalization with nband equal to
                the number of PWs.
        """
        work = cls(manager=manager)
        work.green_nband = green_nband
        gwr_task = GWR_TASK.HDIAGO_FULL if green_nband < 0 else GWR_TASK.HDIAGO
        work.scf_task = work.register_scf_task(scf_input)
        diago_input = scf_input.new_with_vars(optdriver=RUNL.GWR, gwr_task=gwr_task)
        work.diago_task = work.register_gwr_task(diago_input, deps={work.scf_task: "DEN"})

        return work


class _BaseGwrWork(Work):

    @classmethod
    def from_varname_values(cls, varname_values: tuple, gwr_template: AbinitInput,
                            den_node: Node, wfk_node: Node,
                            manager: TaskManager = None):
        """
        Args:
            varname_values:
            gwr_template:
            den_node:
            wfk_node:
            manager:

        Example:

            varname_values = ("nband", [8, 12, 14])

            Work.from_varname_values(varname_values, gwr_template, den_node, wfk_node)

            or:

            varname_values = [
                ("nband", [8, 12]),
                ("gwr_ntau", [6, 8]),
            ]

            Work.from_varname_values(varname_values, gwr_template, den_node, wfk_node)
        """
        work = cls(manager=manager)
        for new_inp in gwr_template.news_varname_values(varname_values):
            work.register_gwr_task(new_inp, deps={den_node: "DEN", wfk_node: "WFK"})

        return work


class GwrSigmaConvWork(_BaseGwrWork):
    """
    This work performs multiple QP calculations with the GWR code
    and produces `xlsx` files in its `outdata` directory
    with the QP results obtained with the different parameters.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GwrSigmaConvWork
    """

    def on_all_ok(self):
        """
        Generate `qp_dirgaps.xlsx` file in the `outdata` directory.
        """
        gwr_files = self.get_all_outdata_files_with_ext("_GWR.nc")
        with GwrRobot.from_files(gwr_files) as gwr_robot:
            dirgaps_df = gwr_robot.get_dirgaps_dataframe()
            dirgaps_df.to_excel(self.outdir.path_in("dirgaps.xlsx"))
            qpdata_df = gwr_robot.get_dataframe()
            qpdata_df.to_excel(self.outdir.path_in("qpdata.xlsx"))

            with gwr_robot.get_pyscript(self.outdir.path_in("gwr_robot.py")) as script:
                script.add_text("""

#robot.plot_selfenergy_conv(spin=0, kpoint=?, band=?, axis="wreal", sortby=None, hue=None)

#robot.plot_qpgaps_convergence(self, qp_kpoints="all", qp_type="qpz0", sortby=None, hue=None)
""")

        return super().on_all_ok()


#class GwrRpaEneConvWork(_BaseGwrWork):
#    """
#    This work computes the RPA energy for different number
#    of points in the minimax mesh.
#    """
#
#    def on_all_ok(self):
#        """
#        """
#        gwr_files = self.get_all_outdata_files_with_ext("_GWR.nc")
#        with GwrRobot.from_files(gwr_files) as robot:
#            df = robot.get_rpaene_dataframe(with_params=True)
#            df.to_excel(self.outdir.path_in("rpaene.xlsx"))
#
#        return super().on_all_ok()
