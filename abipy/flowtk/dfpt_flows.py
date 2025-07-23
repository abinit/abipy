# coding: utf-8
"""
Flow subclasses related to DFPT.
"""
from __future__ import annotations

import dataclasses
import numpy as np

from abipy.tools.typing import TYPE_CHECKING, PathLike
from abipy.tools.serialization import Serializable
from abipy.dfpt.ddb import DdbRobot
from .works import Work, PhononWork
from .flows import Flow

if TYPE_CHECKING:  # needed to avoid circular imports
    from abipy.abio.inputs import AbinitInput


class ConvBecsEpsFlow(Flow):
    """
    This flow performs a convergence study for the Born effective charges,
    the electronic dielectric tensor and the dynamical quadrupoles in semiconductors.
    with respect to the k-point sampling.
    """

    @classmethod
    def from_scf_input(cls,
                       workdir: PathLike,
                       scf_input: AbinitInput,
                       ngkpt_list: list,
                       with_becs: bool = True,
                       with_quad: bool = True,
                       with_flexoe: bool = False,
                       manager=None) -> ConvBecsEpsinfFlow:
        """
        Build a flow for convergence studies wrt ngkpt from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            workdir: Working directory of the flow.
            scf_input: |AbinitInput| for GS-SCF run used as template to generate the other inputs.
            ngkpt_list: List of ngkpt values defining the k-meshes used for DFPT calculations at q = Gamma.
            with_becs: Activate calculation of Electric field and Born effective charges.
            with_quad: Activate calculation of dynamical quadrupoles. Require `with_becs`
                Note that only selected features are compatible with dynamical quadrupoles.
                Please consult <https://docs.abinit.org/topics/longwave/>
            with_flexoe: Activate computation of Flexo-electric tensor.
            manager: |TaskManager| instance. Use default if None.
        """
        flow = cls(workdir=workdir, manager=manager)

        work = Work()
        flow.scf_task = work.register_scf_task(scf_input)
        flow.register_work(work)

        # Save input parameters in flow.
        flow.ngkpt_list = np.reshape(ngkpt_list, (-1, 3))
        flow.with_becs = with_becs
        flow.with_quad = with_quad
        flow.with_flexoe = with_flexoe

        return flow

    def on_all_ok(self):
        """
        This method is called when all the works in the flow have reached S_OK.
        """
        self.on_all_ok_num_calls += 1

        if self.on_all_ok_num_calls == 1:
            # Register NSCF tasks with different k-meshes and only occupied states (TODO)
            self.nscf_tasks = []
            work = self[0]
            for ngkpt in self.ngkpt_list:
                new_input = self.scf_task.input.deepcopy()
                new_input.pop_tolerances()
                new_input.set_vars(ngkpt=ngkpt, iscf=-2, tolwfr=1e-20)
                nscf_task = work.register_nscf_task(new_input, deps={self.scf_task: "DEN"})
                self.nscf_tasks.append(nscf_task)

            self.allocate(build=True)
            return False

        if self.on_all_ok_num_calls == 2:
            # Register DFPT works at Gamma with different k-meshes.
            # Each work reads the WFK file produced by a nscf_task.
            self.conv_works = []
            qpoints = (0, 0, 0)
            for nscf_task in self.nscf_tasks:
                # Important: Remove iscf -2 from the input before passing it to from_scf_task
                nscf_task.input.set_vars(iscf=None, irdden=None)
                work = PhononWork.from_scf_task(nscf_task,
                    qpoints, is_ngqpt=False,
                    with_becs=self.with_becs, with_quad=self.with_quad,
                    with_flexoe=self.with_flexoe, with_dvdb=False,
                    #tolerance=None, ddk_tolerance=None, ndivsm=0, qptopt=1,
                    #prtwf=-1, prepgkk=0, manager=None,
                    )
                self.register_work(work)
                self.conv_works.append(work)

            self.allocate(build=True)
            return False

        if self.on_all_ok_num_calls == 3:
            # Write json file with metadata and the location of the DDB file.
            data = dict(ngkpt_list=self.ngkpt_list,
                        with_becs=self.with_becs,
                        with_quad=self.with_quad,
                        with_flexoe=self.with_flexoe,
            )
            data["ddb_paths"] = [work.outdir.path_in("out_DDB") for work in self.conv_works]

            json_filepath = self.outdir.path_in("ConvBecsEpsResults.json")
            ConvBecsEpsResults(**data).json_write(json_filepath, indent=4)

        return True


@dataclasses.dataclass(kw_only=True)
class ConvBecsEpsResults(Serializable):
    """
    Stores the paths to the DDB files produced with different k-meshes.
    """
    with_becs: bool
    with_quad: bool
    with_flexoe: bool

    ngkpt_list: np.ndarray
    ddb_paths: list[str]  # Paths to DDB files.

    def get_ddb_robot(self) -> DdbRobot:
        """
        Return a DdbRobot instance that can be used to perform convergence studies.
        """
        label_file = {str(ngkpt): ddb_path for ngkpt, ddb_path in zip(self.ngkpt_list, self.ddb_paths, strict=True)}
        ddb_robot = DdbRobot.from_label_file_dict(label_file)
        return ddb_robot
