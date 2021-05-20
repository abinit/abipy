# coding: utf-8
"""Work subclasses related to effective mass calculations."""

import numpy as np
import os

from monty.json import jsanitize
from abipy.core.kpoints import build_segments
from .nodes import Node
from .works import Work, PhononWork
from .flows import Flow


def _get_red_dirs_from_opts(red_dirs, cart_dirs, reciprocal_lattice):
    """
    Helper function to compute the list of directions from user input. Return numpy array.
    """
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
    Useful for cases such as NC+SOC where DFPT is not implemented or if one is interested
    in non-parabolic behaviour.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EffMassLineWork
    """

    @classmethod
    def from_scf_input(cls, scf_input, k0_list, step=0.01, npts=15,
                       red_dirs=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], cart_dirs=None,
                       den_node=None, manager=None):
        """
        Build the Work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            k0_list: List with the reduced coordinates of the k-points where effective masses are wanted.
            step: Step for finite difference in Angstrom^-1
            npts: Number of points sampled around each k-point for each direction.
            red_dirs: List of reduced directions used to generate the segments passing through the k-point
            cart_dirs: List of Cartesian directions used to generate the segments passing through the k-point
            den_node: Path to the DEN file or Task object producing a DEN file.
                Can be used to avoid the initial SCF calculation if a DEN file is already available.
                If None, a GS calculation is performed.
            manager: |TaskManager| instance. Use default if None.
        """
        if npts < 3:
            raise ValueError("Number of points: `%s` should be >= 3 for finite differences" % npts)

        # Define list of directions from user input.
        reciprocal_lattice = scf_input.structure.lattice.reciprocal_lattice
        all_red_dirs = _get_red_dirs_from_opts(red_dirs, cart_dirs, reciprocal_lattice)
        kpts = build_segments(k0_list, npts, step, all_red_dirs, reciprocal_lattice)

        # Now build NSCF input with explicit list of k-points.
        nscf_input = scf_input.make_nscf_kptopt0_input(kpts)

        new = cls(manager=manager)

        # Need to perform SCF run if DEN file is not available.
        scf_task = new.register_scf_task(scf_input) if den_node is None else Node.as_node(den_node)

        new.register_nscf_task(nscf_input, deps={scf_task: "DEN"})
        return new


class EffMassDFPTWork(Work):
    """
    Work for the computation of effective masses with DFPT.
    Requires explicit list of k-points and range of bands.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EffMassDFPTWork
    """

    @classmethod
    def from_scf_input(cls, scf_input, k0_list, effmass_bands_f90, ngfft=None, den_node=None, manager=None):
        """
        Build the Work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            k0_list: List with the reduced coordinates of the k-points where effective masses are wanted.
            effmass_bands_f90: (nkpt, 2) array with band range for effmas computation.
                WARNING: Assumes Fortran convention with indices starting from 1.
            ngfft: FFT divisions (3 integers). Used to enforce the same FFT mesh in the NSCF run as the one used for GS.
            den_node: Path to the DEN file or Task object producing a DEN file.
                Can be used to avoid the initial SCF calculation if a DEN file is already available.
                If None, a GS calculation is performed.
            manager: |TaskManager| instance. Use default if None.
        """
        multi = scf_input.make_dfpt_effmass_inputs(k0_list, effmass_bands_f90, ngfft=ngfft)
        nscf_input, effmass_input = multi[0], multi[1]

        new = cls(manager=manager)

        # Important: keep a reference to the Frohlich input that can be used to run EPH calculations if needed.
        new.frohlich_input = multi[2]

        # Need to perform SCF run if DEN file is not available
        scf_task = new.register_scf_task(scf_input) if den_node is None else Node.as_node(den_node)

        nscf_task = new.register_nscf_task(nscf_input, deps={scf_task: "DEN"})
        new.register_effmass_task(effmass_input, deps={scf_task: "DEN", nscf_task: "WFK"})
        return new


class EffMassAutoDFPTWork(Work):
    """
    Work for the automatic computation of effective masses with DFPT.
    Band extrema are automatically detected by performing a NSCF calculation
    along a high-symmetry k-path with ndivsm.
    Requires more computation that EffMassWork since input variables (kpoints and band range)
    are computed at runtime.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: EffMassAutoDFPTWork
    """

    @classmethod
    def from_scf_input(cls, scf_input, ndivsm=15, tolwfr=1e-20, den_node=None, manager=None):
        """
        Build the Work from an |AbinitInput| representing a GS-SCF calculation.

        Args:
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            ndivsm: if > 0, it's the number of divisions for the smallest segment of the path (Abinit variable).
                if < 0, it's interpreted as the pymatgen `line_density` parameter in which the number of points
                in the segment is proportional to its length. Typical value: -20.
                This option is the recommended one if the k-path contains two high symmetry k-points that are very close
                as ndivsm > 0 may produce a very large number of wavevectors.
            tolwfr: Tolerance on residuals for NSCF calculation
            den_node: Path to the DEN file or Task object producing a DEN file.
                Can be used to avoid the initial SCF calculation if a DEN file is already available.
                If None, a GS calculation is performed.
            manager: |TaskManager| instance. Use default if None.
        """
        if scf_input.get("nsppol", 1) != 1:
            raise NotImplementedError("Magnetic semiconductors with nsppol = 2 are not implemented!")

        new = cls(manager=manager)
        # Keep a copy of the initial input.
        new.scf_input = scf_input.deepcopy()

        # Need SCF run to get DEN file.
        if den_node is None:
            new.den_node = new.register_scf_task(new.scf_input)
        else:
            new.den_node = Node.as_node(den_node)

        # Perform NSCF run along k-path that will be used to find band extrema.
        bands_input = scf_input.make_ebands_input(ndivsm=ndivsm, tolwfr=tolwfr)
        new.bands_task = new.register_nscf_task(bands_input, deps={new.den_node: "DEN"})

        return new

    def on_all_ok(self):
        """
        This method is called once the `Work` is completed i.e. when all tasks have reached status S_OK.
        Here, we read the band structure from GSR to find the position of the band edges and use these values
        to generate a new work for effective masses with DFPT.
        """
        with self.bands_task.open_gsr() as gsr:
            ebands = gsr.ebands
            # Warning: Assuming semiconductor with spin-unpolarized band energies.
            # At present, Abinit input variables do not take into account nsppol.
            ebands.set_fermie_to_vbm()
            # Find k0_list and effmass_bands_f90
            k0_list, effmass_bands_f90 = ebands.get_kpoints_and_band_range_for_edges()
            den_ngfft = gsr.reader.read_ngfft3()

        # Create the work for effective mass computation with DFPT and add it to the flow.
        # Keep a reference in generated_effmass_dfpt_work.
        work = EffMassDFPTWork.from_scf_input(self.scf_input, k0_list, effmass_bands_f90,
                                              ngfft=den_ngfft, den_node=self.den_node)

        self.generated_effmass_dfpt_work = work
        self.flow.register_work(work)
        self.flow.allocate()
        self.flow.build_and_pickle_dump()
        self.flow.finalized = False

        return super().on_all_ok()


class FrohlichZPRFlow(Flow):
    """
    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: FrohlichZPRFlow
    """

    @classmethod
    def from_scf_input(cls, workdir, scf_input, ddb_node=None, ndivsm=15, tolwfr=1e-20, metadata=None, manager=None):
        """
        Build the Work from an |AbinitInput| representing a GS-SCF calculation.
        Final results are stored in the "zprfrohl_results.json" in the outdata directory of the flow.

        Args:
            workdir: Working directory.
            scf_input: |AbinitInput| for GS-SCF used as template to generate the other inputs.
            ddb_node: Path to an external DDB file that is used to avoid the calculation of BECS/eps_inf and phonons.
                If None, a DFPT calculation is automatically performed by the flow.
            ndivsm: Number of divisions used to sample the smallest segment of the k-path.
            tolwfr: Tolerance on residuals for NSCF calculation
            manager: |TaskManager| instance. Use default if None.
            metadata: Dictionary with metadata addeded to the final JSON file.
        """
        new = cls(workdir=workdir, manager=manager)
        new.metadata = jsanitize(metadata) if metadata is not None else None

        # Build work for the automatic computation of effective masses.
        new.effmass_auto_work = EffMassAutoDFPTWork.from_scf_input(scf_input, ndivsm=ndivsm, tolwfr=tolwfr)
        new.register_work(new.effmass_auto_work)
        new.scf_task, new.ebands_kpath_task = new.effmass_auto_work[0], new.effmass_auto_work[1]

        if ddb_node is not None:
            new.ddb_node = Node.as_node(ddb_node)
            new.ddb_file_path_if_ddb_node = os.path.abspath(ddb_node)
        else:
            # Compute DDB with BECS and eps_inf.
            new.ddb_file_path_if_ddb_node = None
            ph_work = PhononWork.from_scf_task(new.scf_task, qpoints=[0, 0, 0],
                                               is_ngqpt=False, tolerance=None, with_becs=True,
                                               ddk_tolerance=None)
            new.register_work(ph_work)
            new.ddb_node = ph_work

        return new

    def on_all_ok(self):
        """
        This method is called when all the works in the flow have reached S_OK.
        This method shall return True if the calculation is completed or
        False if the execution should continue due to side-effects such as adding a new work to the flow.
        """
        if self.on_all_ok_num_calls > 0: return True
        self.on_all_ok_num_calls += 1

        work = Work()
        inp = self.effmass_auto_work.generated_effmass_dfpt_work.frohlich_input
        wfk_task = self.effmass_auto_work.generated_effmass_dfpt_work[0]
        self.effmass_task = self.effmass_auto_work.generated_effmass_dfpt_work[1]

        deps = {wfk_task: "WFK", self.ddb_node: "DDB", self.effmass_task: "EFMAS.nc"}
        self.frohl_task = work.register_eph_task(inp, deps=deps)
        self.register_work(work)

        self.allocate()
        self.build_and_pickle_dump()
        return False

    def finalize(self):
        """
        This method is called when the flow is completed.
        Here we write the final results in the "zprfrohl_results.json" file
        in the `outdata` directory of the flow
        """
        d = {}
        if self.metadata is not None: d.update({"metadata": self.metadata})

        # Add GS results.
        with self.scf_task.open_gsr() as gsr:
            d["gsr_scf_path"] = gsr.filepath
            d["pressure_GPa"] = float(gsr.pressure)
            d["max_force_eV_Ang"] = float(gsr.max_force)
            d["structure"] = gsr.structure

        # Add NSCF band structure.
        with self.ebands_kpath_task.open_gsr() as gsr:
            d["gsr_nscf_kpath"] = gsr.filepath
            gsr.ebands.set_fermie_to_vbm()
            d["ebands_kpath"] = gsr.ebands
            #d["ebands_info"] = gsr.ebands.get_dict4pandas(with_geo=False, with_spglib=False)

        # TODO
        # Extract results from run.abo
        #d["effmass_results"] = self.effmass_task.yaml_parse_results()
        #d["frohl_results"] = self.frohl_task.yaml_parse_results()

        # Add epsinf, e0, BECS, alpha and DDB as string.
        from abipy import abilab
        if self.ddb_file_path_if_ddb_node is not None:
            ddb_filepath = self.ddb_file_path_if_ddb_node
        else:
            ddb_filepath = self.ddb_node.outdir.path_in("out_DDB")

        with abilab.abiopen(ddb_filepath) as ddb:
            d["ddb_path"] = ddb.filepath
            d["ddb_string"] = ddb.get_string()
            epsinf, becs = ddb.anaget_epsinf_and_becs(chneut=1)
            gen = ddb.anaget_dielectric_tensor_generator(asr=2, chneut=1, dipdip=1)
            eps0 = gen.tensor_at_frequency(w=0.0)
            d["becs"] = becs
            d["epsinf_cart"] = epsinf
            # FIXME Complex is not supported by JSON.
            d["eps0_cart"] = eps0.real

        #print(d)
        abilab.mjson_write(d, self.outdir.path_in("zprfrohl_results.json"), indent=4)

        return super().finalize()
