# coding: utf-8
"""Work subclasses related to effective mass calculations."""

import numpy as np
import os

from abipy.core.kpoints import Kpoint
from .nodes import Node
from .works import Work, PhononWork
from .flows import Flow


def _get_red_dirs_from_opts(red_dirs, cart_dirs, reciprocal_lattice):
    """Helper function to compute list of directions from user input. Return numpy array."""
    all_red_dirs = []
    if red_dirs is not None:
        all_red_dirs.extend(np.reshape(red_dirs, (-1, 3)))

    if cart_dirs is not None:
        frac_coords = reciprocal_lattice.get_fractional_coords(np.reshape(cart_dirs, (-1, 3)))
        all_red_dirs.extend(frac_coords)

    return np.reshape(all_red_dirs, (-1, 3))


class EffMassLineWork(Work):
    """
    Work for the computation of effective masses via finite differences along a k-line.
    Useful for cases such as NC+SOC where DFPT is not coded or for debugging purposes.
    """

    @classmethod
    def from_scf_input(cls, scf_input, k0_list, step=0.01, npts=15,
                       red_dirs=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], cart_dirs=None,
                       den_node=None, manager=None):
        """
        Build the Work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            k0_list: List with the reduced coordinates of the k-point where effective masses are wanted.
            step: Step in Angstrom^-1
            npts: Number of points sampled around each kpoint for each direction.
            red_dirs: List of reduced directions used to generate the segments passing through the k-point
            cart_dirs: List of Cartesian directions used to generate the segments passing through the k-point
            den_node: Path to the DEN file or Task object producing a DEN file.
                Can be used to avoid the initial SCF calculation if a DEN file is already available.
            manager: |TaskManager| instance. Use default if None.
        """
        if npts < 3:
            raise ValueError("Number of points: `%s` should be >= 3 for finite differences" % npts)

        reciprocal_lattice = scf_input.structure.lattice.reciprocal_lattice

        # Define list of directions from user input.
        all_red_dirs = _get_red_dirs_from_opts(red_dirs, cart_dirs, reciprocal_lattice)

        kpts = []
        for kpoint in np.reshape(k0_list, (-1, 3)):
            kpoint = Kpoint.as_kpoint(kpoint, reciprocal_lattice)
            # Build segment passing through this kpoint (work in Cartesian coords)
            for rdir in all_red_dirs:
                bvers = reciprocal_lattice.matrix.T @ rdir
                bvers /= np.sqrt(np.dot(bvers, bvers))
                k0 = kpoint.cart_coords - bvers * (npts // 2) * step
                for ii in range(npts):
                    kc = k0 + ii * bvers * step
                    kpts.append(kc)

        # Cart --> Frac then build NSCF input with explicit list of k-points.
        kpts = reciprocal_lattice.get_fractional_coords(np.reshape(kpts, (-1, 3)))
        nscf_input = scf_input.make_nscf_kptopt0_input(kpts)

        new = cls(manager=manager)

        # Need to perform SCF run if DEN file is not available.
        scf_task = new.register_scf_task(scf_input) if den_node is None else Node.as_node(den_node)

        new.register_nscf_task(nscf_input, deps={scf_task: "DEN"})
        return new


class EffMassDFPTWork(Work):
    """
    Work for the computation of effective masses with DFPT.
    Requires explicit list of k-points and band range.
    """

    @classmethod
    def from_scf_input(cls, scf_input, k0_list, effmass_bands_f90,
                       red_dirs=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], cart_dirs=None,
                       den_node=None, manager=None):
        """
        Build the Work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            k0_list: List with the reduced coordinates of the k-point where effective masses are wanted.
            effmass_bands_f90: (nkpt, 2) array with band range for effmas computation.
                WARNING: Uses Fortran convention starting from 1.
            red_dirs: List of reduced directions used to generate the segments passing through the k-point
            cart_dirs: List of Cartesian directions used to generate the segments passing through the k-point
            den_node: Path to the DEN file or Task object producing a DEN file.
                Can be used to avoid the initial SCF calculation if a DEN file is already available.
            manager: |TaskManager| instance. Use default if None.
        """
        reciprocal_lattice = scf_input.structure.lattice.reciprocal_lattice

        # Define list of directions from user input.
        all_red_dirs = _get_red_dirs_from_opts(red_dirs, cart_dirs, reciprocal_lattice)

        multi = scf_input.make_dfpt_effmass_input(k0_list, effmass_bands_f90)
        nscf_input, effmass_input = multi[0], multi[1]

        new = cls(manager=manager)

        # Important: keep a reference to the Frohlich input that can be used to run EPH calculations if needed.
        new.frohlich_input = multi[2]

        # Need to perform SCF run if DEN file is not available
        scf_task = new.register_scf_task(scf_input) if den_node is None else Node.as_node(den_node)

        nscf_task = new.register_nscf_task(nscf_input, deps={scf_task: "DEN"})
        new.register_effmass_task(effmass_input, deps={nscf_task: ["DEN", "WFK"]})
        return new


class EffMassAutoDFPTWork(Work):
    """
    Work for the automatic computation of effective masses with DFPT.
    Band extrema are automatically detected by performing a NSCF calculation along a high-symmetry k-path.
    Require more computation that EffMassWork but input parameters since input variables
    are automatically defined at runtime.
    """

    @classmethod
    def from_scf_input(cls, scf_input, ndivsm=15, tolwfr=1e-20,
                       #red_dirs=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], cart_dirs=None,
                       manager=None):
        """
        Build the Work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            red_dirs: List of reduced directions used to generate the segments passing through the k-point
            cart_dirs: List of Cartesian directions used to generate the segments passing through the k-point
            manager: |TaskManager| instance. Use default if None.
        """
        if scf_input.get("nsppol", 1) != 1:
            raise NotImplementedError("Magnetic semiconductors with nsppol = 2 are not implemented!")

        new = cls(manager=manager)
        # Keep a copy of the initial input.
        new.scf_input = scf_input.deepcopy()

        # Need SCF run to get DEN file.
        new.scf_task = new.register_scf_task(new.scf_input)

        # Perform NSCF run along k-path that will be used to find band extrema.
        bands_input = scf_input.make_bands_input(ndivsm=ndivsm, tolwfr=tolwfr)
        new.bands_task = new.register_nscf_task(bands_input, deps={new.scf_task: "DEN"})

        return new

    def on_all_ok(self):
        """
        This method is called once the `Work` is completed i.e. when all tasks have reached status S_OK.
        Here, we read the band structure from GSR to find the position of the band edges and use these values
        to generate a new work for effective masses with DFPT.
        """
        #print("in on all_ok")
        with self.bands_task.open_gsr() as gsr:
            ebands = gsr.ebands
            # Warning: Assuming semiconductor with spin-unpolarized band energies.
            # At present, Abinit input variables do not take into account nsppol.
            ebands.set_fermie_to_vbm()
            # Find k0_list and effmass_bands_f90
            k0_list, effmass_bands_f90 = ebands.get_kpoints_and_band_range_for_edges()

        # Create the work for effective mass computation with DFPT and add it to the flow.
        # Keep also a reference in generated_effmass_dfpt_work.
        work = EffMassDFPTWork.from_scf_input(self.scf_input, k0_list, effmass_bands_f90,
                                              #red_dirs=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], cart_dirs=None,
                                              den_node=self.scf_task)

        self.generated_effmass_dfpt_work = work
        self.flow.register_work(work)
        self.flow.allocate()
        self.flow.build_and_pickle_dump()
        self.flow.finalized = False

        return super().on_all_ok()


class FrohlichZPRFlow(Flow):
    """
    """

    @classmethod
    def from_scf_input(cls, scf_input, ddb_node=None, ndivsm=15, tolwfr=1e-20,
                       workdir=None, manager=None):
        """
        Build the Work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            ddb_node: Path to an external DDB file that is used to avoid the calculation of BECS/eps_inf and phonons.
                If None, a DFPT calculation is automatically perfored by the flow.
            ndivsm:
            tolwfr:
            workdir:
            manager: |TaskManager| instance. Use default if None.
        """
        new = cls(workdir=workdir, manager=manager)

        # Build work for the automatic computation of effective masses.
        new.effmass_auto_work = EffMassAutoDFPTWork.from_scf_input(scf_input, ndivsm=ndivsm, tolwfr=tolwfr)
        new.register_work(new.effmass_auto_work)
        scf_task = new.effmass_auto_work[0]

        if ddb_node is not None:
            new.ddb_node = Node.as_node(ddb_node)
        else:
            # Compute DDB with BECS and eps_inf.
            becs_work = PhononWork.from_scf_task(scf_task, qpoints=[0, 0, 0],
                                                 is_ngqpt=False, tolerance=None, with_becs=True,
                                                 ddk_tolerance=None)
            new.register_work(becs_work)
            new.ddb_node = becs_work

        new.on_all_ok_num_calls = 0

        return new

    def on_all_ok(self):
        self.on_all_ok_num_calls += 1
        if self.on_all_ok_num_calls > 1:
            print("flow Returning True")
            return True

        work = Work()
        inp = self.effmass_auto_work.generated_effmass_dfpt_work.frohlich_input
        wfk_task = self.effmass_auto_work.generated_effmass_dfpt_work[0]
        effmass_task = self.effmass_auto_work.generated_effmass_dfpt_work[1]
        t = work.register_eph_task(inp, deps={wfk_task: "WFK", self.ddb_node: "DDB", effmass_task: "EFMAS.nc"})

        self.register_work(work)
        self.allocate()
        self.build_and_pickle_dump()
        self.finalized = False
        print("flow Returning False")
        return False
