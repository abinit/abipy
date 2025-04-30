# coding: utf-8
"""
Works and Flows for GWR calculations (GW in supercells).

NB: An Abinit build with Scalapack is required to run GWR.
"""
from __future__ import annotations

import os

from abipy.abio.inputs import AbinitInput, RUNL, GWR_TASK
from abipy.electrons.gwr import GwrRobot
from abipy.tools.iotools import make_executable
from .nodes import Node
from .tasks import TaskManager
from .works import Work


class DirectDiagoWork(Work):
    """
    This work performs the direct diagonalization of the KS Hamiltonian
    using the density produced by a GS run and produces a WFK file with
    empty states in the outdir of the second task.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: DirectDiagoWork
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       green_nband: int,
                       manager: TaskManager=None) -> DirectDiagoWork:
        """
        Build object from an input representing a GS-SCF calculation.

        Args:
            scf_input: Input for the GS-SCF calculation.
            green_nband: Number of bands to compute in the direct diagonalization.
                A negative value activate the full diagonalization with nband
                equal to the number of PWs.
        """
        work = cls(manager=manager)
        work.green_nband = green_nband
        gwr_task = GWR_TASK.HDIAGO_FULL if green_nband < 0 else GWR_TASK.HDIAGO
        work.scf_task = work.register_scf_task(scf_input)
        diago_input = scf_input.new_with_vars(optdriver=RUNL.GWR, gwr_task=gwr_task)
        work.diago_task = work.register_gwr_task(diago_input, deps={work.scf_task: "DEN"})

        return work


class _BaseGWRWork(Work):
    """Base class for GWR works."""

    @classmethod
    def from_varname_values(cls,
                            varname_values: tuple,
                            gwr_template: AbinitInput,
                            den_node: Node,
                            wfk_node: Node,
                            manager: TaskManager = None):
        """
        Generate the work by changing the values of selected variables in a template for GWR calculations.

        Args:
            varname_values: Tuple or List of tuples. See examples below.
            gwr_template: GWR input to be used as template to generate other inputs.
            den_node: The Node who produced the DEN file.
            wfk_node: The Node who produced the WFK file.
            manager: Abipy Task Manager.

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


class GWRSigmaConvWork(_BaseGWRWork):
    """
    This work performs multiple GW calculations with the GWR code
    and produces `xlsx` files in its `outdata` directory
    with the QP results obtained with the different parameters.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GWRSigmaConvWork
    """

    def on_all_ok(self):
        """
        Generate `qp_dirgaps.csv` file in the `outdata` directory.
        """
        with GwrRobot.from_files(self.get_all_outdata_files_with_ext("_GWR.nc")) as gwr_robot:
            dirgaps_df = gwr_robot.get_dirgaps_dataframe()
            dirgaps_df.to_csv(self.outdir.path_in("dirgaps.csv"))
            qpdata_df = gwr_robot.get_dataframe()
            qpdata_df.to_csv(self.outdir.path_in("qpdata.csv"))

            with gwr_robot.get_pyscript(self.outdir.path_in("gwr_robot.py")) as script:
                script.add_text("""
#robot.plot_selfenergy_conv(spin=0, kpoint=?, band=?, axis="wreal", sortby=None, hue=None)
#robot.plot_qpgaps_convergence(qp_kpoints="all", qp_type="qpz0", sortby=None, hue=None)
""")

        return super().on_all_ok()



class GWRChiCompareWork(_BaseGWRWork):
    """
    This work computes the irreducibile polarizability along the imaginary axis
    using the GWR code and the quartic-scaling algorithm using the same minimax mesh
    so that one can compare the two quantities.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: GWRChiCompareConvWork
    """

    @classmethod
    def from_scf_input(cls,
                       scf_input: AbinitInput,
                       gwr_ntau: int,
                       nband: int,
                       ecuteps: float,
                       den_node: Node,
                       wfk_node: Node,
                       gwr_kwargs: dict | None = None,
                       scr_kwargs: dict | None = None,
                       manager: TaskManager = None):
        """
        Build Work from an input for GS-SCF calculation

        Args:
            scf_input: input for GS run.
            gwr_ntau: Number of points in minimax mesh.
            nband: Number of bands to build G and chi.
            ecuteps: Cutoff energy for chi.
            den_node: The Node who produced the DEN file.
            wfk_node: The Node who produced the WFK file.
            gwr_kwargs: Extra kwargs used to build the GWR input.
            scr_kwargs: Extra kwargs used to build the SCR input.
            manager: Abipy Task Manager.
        """
        gwr_input = scf_input.make_gwr_qprange_input(gwr_ntau=gwr_ntau, nband=nband,
                                                     ecuteps=ecuteps, gwr_task=GWR_TASK.CHI0)
        gwr_input.set_vars(prtsuscep=1, iomode=3)
        if gwr_kwargs is not None: gwr_input.set_vars(**gwr_kwargs)

        chi_input = scf_input.new_with_vars(optdriver=3,
                                            gwcalctyp=1, # Analytic continuation.
                                            nfreqim=gwr_ntau,
                                            ecuteps=ecuteps,
                                            nband=nband,
                                            prtsuscep=1,
                                            iomode=3,
                                            userie=4242, # Magic number to use the minimax mesh in the SCR driver.
                                                         # with nfreqim points.
                                            )
        if scr_kwargs is not None: chi_input.set_vars(**scr_kwargs)

        work = cls(manager=manager)
        work.register_scr_task(chi_input, deps={wfk_node: "WFK"})
        work.register_gwr_task(gwr_input, deps={den_node: "DEN", wfk_node: "WFK"})
        return work

    def on_all_ok(self):
        """
        Write python script to compare chi0 matrix elements.
        """
        sus_path = self[0].outdir.path_in("out_SUS.nc")
        tchi_path = self[1].outdir.path_in("out_TCHIM.nc")

        py_text = f"""\
#!/usr/bin/env python

import seaborn as sns
sns.set(context="paper", style='darkgrid', palette='deep',
       font='sans-serif', font_scale=0.8, color_codes=False, rc=None)

from abipy.electrons.gwr import TchimVsSus
o = TchimVsSus("{tchi_path}", "{sus_path}")

gpairs = [
    ((0, 0, 0), (0, 0, 0)),
    ((0, 0, 0), (1, 0, 0)),
    ((1, 0, 0), (1, 0, 0)),
    ((1, 0, 0), (0, 1, 0)),
    ((0, 0, 0), (2, 0, 0)),
]

qpoint_list = [
    [0, 0, 0],
    [0.5, 0, 0],
]

o.expose_qpoints_gpairs(qpoint_list, gpairs, exposer="mpl")

#o.plot_matdiff(qpoint=0, iw_index=0)
"""
        py_path = os.path.join(self.workdir, "plot_tchi_sus.py")
        with open(py_path, "wt") as fh:
            fh.write(py_text)
        make_executable(py_path)

        return super().on_all_ok()


class GWRRPAConvWork(_BaseGWRWork):
    """
    This work computes the RPA correlated energy for different number
    of points in the minimax mesh.
    """

    @classmethod
    def from_scf_input_ntaus(cls, scf_input: AbinitInput,
                             gwr_ntau_list: list,
                             nband: int,
                             ecuteps: float,
                             den_node: Node,
                             wfk_node: Node,
                             gwr_kwargs: dict | None = None,
                             manager: TaskManager = None):
        """
        Build Work from an input for GS-SCF calculation

        Args:
            scf_input: input for GS run.
            gwr_ntau_list: List with number of points in minimax mesh.
            nband: Number of bands to build G and chi.
            ecuteps: Cutoff energy for chi.
            den_node: The Node who produced the DEN file.
            wfk_node: The Node who produced the WFK file.
            gwr_kwargs: Extra kwargs used to build the GWR input.
            manager: Abipy Task Manager.
        """
        work = cls(manager=manager)
        for gwr_ntau in gwr_ntau_list:
            gwr_input = scf_input.make_gwr_qprange_input(gwr_ntau=gwr_ntau, nband=nband,
                                                         ecuteps=ecuteps, gwr_task=GWR_TASK.RPA_ENERGY)
            if gwr_kwargs is not None: gwr_input.set_vars(**gwr_kwargs)
            work.register_gwr_task(gwr_input, deps={den_node: "DEN", wfk_node: "WFK"})

        return work

    def on_all_ok(self):
        """
        Generate `gwr_robot.py` script in the `outdata` directory.
        """
        with GwrRobot.from_files(self.get_all_outdata_files_with_ext("_GWR.nc")) as gwr_robot:
            with gwr_robot.get_pyscript(self.outdir.path_in("gwr_robot.py")) as script:
                script.add_text("""
df = robot.get_rpa_ene_dataframe()
print(df)
""")

        return super().on_all_ok()
