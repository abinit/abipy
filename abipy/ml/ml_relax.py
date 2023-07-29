"""
"""
from __future__ import annotations

import os
import time
import tempfile
import json
import numpy as np

from pathlib import Path
from typing import Any
#from monty.string import marquee, list_strings # is_string,
from monty.json import MontyEncoder
from ase.atoms import Atoms
from ase.optimize import BFGS
from ase.calculators.abinit import Abinit
from ase.constraints import ExpCellFilter
from ase.stress import full_3x3_to_voigt_6_stress, voigt_6_to_full_3x3_stress
from abipy.core.abinit_units import eV_Ha, Ang_Bohr
from abipy.tools.iotools import workdir_with_prefix
from abipy.ml.aseml import print_atoms, get_atoms, get_structure, CalcBuilder, scompare_two_atoms
from abipy.dynamics.hist import HistFile
from abipy.flowtk import PseudoTable


class RelaxationProfiler:

    def __init__(self, atoms: Any, pseudos, relax_mode: str, fmax: float, verbose: int):
        """
        """
        self.atoms = get_atoms(atoms)
        self.relax_mode = relax_mode
        assert self.relax_mode in ("ions", "cell")
        self.fmax = fmax
        self.verbose = verbose

        pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(get_structure(self.atoms))
        pp_paths = [p.filepath for p in pseudos]

        hints = [p.hint_for_accuracy("normal") for p in pseudos]
        ecut = max(h.ecut for h in hints) * 27.3  # In ASE this is in eV (don't know why!)
        #pawecutdg = max(h.pawecutdg for h in hints) if pseudos.allpaw else None

        self.gs_kwargs = dict(
            ecut=ecut,
            tolvrs=1e-8,
            kpts=[2, 2, 2],
            nstep=100,
            expert_user=1,   # Ignore warnings (chksymbreak, chksymtnons, chkdilatmx)
            autoparal=1,
            paral_kgb=1,
            pseudos=pp_paths,
        )

    def run(self, nn_name="m3gnet", workdir=None, prefix=None):
        """
        """
        workdir = workdir_with_prefix(workdir, prefix)
        print(f"Running in {str(workdir)} ...")

        # Run fully ab-initio relaxation with abinit.
        # TODO: Fix issue with ixc set by ASE.
        relax_kwargs = dict(
            ecutsm=0.5,     # Smoothing PW cutoff energy (mandatory for cell optimization)
            #ionmov=2,
            ionmov=22,
            #ionmov=28,     # activate i-pi/socket mode
            optcell=0 if self.relax_mode == "ions" else 2,
            tolmxf=self.fmax * eV_Ha * Ang_Bohr,
            ntime=100,
        )
        relax_kwargs.update(**self.gs_kwargs)

        print("Begin relaxation with ABINIT and tolmxf:", relax_kwargs["tolmxf"])
        abinit_relax = Abinit(directory=workdir / "abinit_relax", **relax_kwargs)
        abi_atoms = self.atoms.copy()
        t_start = time.time()
        abi_forces = abinit_relax.get_forces(atoms=abi_atoms)
        abi_fmax = np.sqrt((abi_forces ** 2).sum(axis=1).max())
        abi_converged = abi_fmax < self.fmax
        if self.verbose:
            print("abi_fmax:", abi_fmax, ", fmax:", self.fmax, ", converged:", abi_converged)
            #print("abi_forces:\n", abi_forces)
        t_end = time.time()

        with HistFile(abinit_relax.directory / "abinito_HIST.nc") as hist:
            abi_nsteps = hist.num_steps
            #print(hist)
        print('ABINIT relaxation completed in %2.4f sec after nsteps %d\n' % (t_end - t_start, abi_nsteps))

        # Run hybrid relaxation (ML + abinit)
        if self.verbose:
            ml_calc = CalcBuilder(nn_name).get_calculator()
        else:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ml_calc = CalcBuilder(nn_name).get_calculator()

        count, abiml_nsteps, ml_nsteps = 0, 0, 0
        t_start = time.time()
        while count < 10:
            count += 1
            # Compute ab-initio forces and check for convergence.
            abinit = Abinit(directory=workdir / "abiml_relax", **self.gs_kwargs)
            # This one seems to be needed to get updated forces but don't know why!!
            abinit.use_cache = False
            abiml_nsteps += 1
            abi_forces = abinit.get_forces(atoms=self.atoms)
            abi_stress = abinit.get_stress(atoms=self.atoms)
            abi_stress = voigt_6_to_full_3x3_stress(abi_stress)
            abi_fmax = np.sqrt((abi_forces ** 2).sum(axis=1).max())
            abi_converged = abi_fmax < self.fmax
            #print("Iteration count:", count, "abi_fmax:", abi_fmax, ", fmax:", self.fmax, ", converged:", abi_converged)
            #print("abi_stress", abi_stress)
            #print_atoms(atoms, cart_forces=abi_forces)
            if abi_converged:
                print("Converged with ABINIT forces. breaking now!")
                break

            #if abi_fmax > 100 * self.fmax:
            #if (abi_fmax > 100 * self.fmax) and False:
            if True:
                # Compute ML forces and set delta forces in the ML calculator.
                ml_calc.set_delta_forces(None)
                ml_forces = ml_calc.get_forces(atoms=self.atoms)
                delta_forces = abi_forces - ml_forces
                #delta_forces = None
                ml_calc.set_delta_forces(delta_forces)
                if self.relax_mode == "cell":
                    ml_calc.set_delta_stress(None)
                    ml_stress = ml_calc.get_stress(atoms=self.atoms)
                    delta_stress = abi_stress - ml_stress
                    ml_calc.set_delta_stress(delta_stress)

                self.atoms.calc = ml_calc
                has_ml_calc = True
            else:
                self.atoms.calc = abinit
                has_ml_calc = False

            #print("Begin optimization with has_ml_calc:", has_ml_calc)
            opt_kws = dict(
                trajectory=str(abinit.directory / f"opt_{count}.traj"),
                logfile=str(abinit.directory / f"log_{count}"),
            )

            if self.relax_mode == "cell":
                opt = BFGS(ExpCellFilter(self.atoms), **opt_kws)
            else:
                opt = BFGS(self.atoms, **opt_kws)

            opt.run(fmax=self.fmax, steps=None)
            if has_ml_calc:
                ml_nsteps += opt.nsteps
            else:
                abiml_nsteps += opt.nsteps

            #print(scompare_two_atoms("ABINIT relaxed", abi_atoms, "ABIML", self.atoms))
            opt_converged = opt.converged()
            if opt_converged and not has_ml_calc:
                print("ASE optimizer converged with ABINIT forces. breaking now!")
                break

        t_end = time.time()
        print_atoms(self.atoms, title="Atoms after relaxation:")
        print(f"{ml_nsteps=}, {abiml_nsteps=}")
        print('ML + ABINIT Relaxation completed in %2.4f sec\n' % (t_end - t_start))

        # Write json file with output results.
        with open(workdir / "data.json", "wt") as fh:
            data = dict(
                abi_nsteps=abi_nsteps,
                abiml_nsteps=abiml_nsteps,
                ml_nsteps=ml_nsteps,
                abi_relaxed_structure=get_structure(abi_atoms),
                abiml_relaxed_structure=get_structure(self.atoms),
            )
            json.dump(data, fh, indent=4, cls=MontyEncoder)
            # Print JSON data
            #print("")
            #print(marquee("Results", mark="="))
            #print(json.dumps(data, indent=4, cls=MontyEncoder), end="\n")


if __name__ == "__main__":
    from abipy.flowtk.psrepos import get_repo_from_name
    pseudos = get_repo_from_name("ONCVPSP-PBE-SR-PDv0.4").get_pseudos("standard")
    from ase.build import bulk
    atoms = bulk('Si')
    atoms.rattle(stdev=0.1, seed=42)
    prof = RelaxationProfiler(atoms, pseudos, relax_mode="ions", fmax=0.001, verbose=0)
    prof.run()
