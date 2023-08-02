#!/usr/bin/env python
"""
Script to perform structural relaxations with ML + ABINIT
"""
from __future__ import annotations


import sys
import os
import click
import numpy as np

from ase.atoms import Atoms
from ase.build import bulk
from abipy.core.structure import Structure
from abipy.flowtk.psrepos import get_repo_from_name
from abipy.ml.aseml import get_atoms, ase_optimizer_cls
from abipy.ml.ml_relax import RelaxationProfiler


ASE_OPTIMIZERS = ase_optimizer_cls("__all__")


def add_relax_opts(f):
    """Add CLI options for structural relaxations."""
    # fmax (float): total force tolerance for relaxation convergence.
    # Here fmax is a sum of force and stress forces. Defaults to 0.1.
    f = click.option("--relax-mode", "-r", default="cell", show_default=True, type=click.Choice(["no", "ions", "cell"]),
                     help='Relaxation mode.')(f)
    f = click.option("--fmax", default=0.005, type=float, show_default=True,
                     help='Stopping criterion in eV/A')(f)
    f = click.option("--steps", default=500, type=int, show_default=True,
                     help="Max number of steps for ASE relaxation.")(f)
    f = click.option("--optimizer", "-o", default="BFGS", show_default=True, type=click.Choice(ASE_OPTIMIZERS),
                     help="ASE optimizer class.")(f)
    return f


def add_workdir_verbose_opts(f):
    """Add workdir and verbose options to CLI subcommand."""
    f = click.option("--workdir", "-w", default=None, type=str,
                      help="Working directory. If not specified, a temporary directory is created.")(f)
    f = click.option('-v', '--verbose', count=True, help="Verbosity level")(f)
    return f


@click.command()
@click.argument("filepath", type=str)
@click.option("--nn-name", "-nn", default="m3gnet", show_default=True,
              type=click.Choice(["m3gnet_old", "m3gnet", "chgnet"]), help='ML potential')
@add_relax_opts
@click.option("-k", "--kppa", default=1000, type=float, show_default=True,
              help="Defines the sampling of BZ mesh (k-points per atom)")
@click.option("--rattle", default=0.0, type=float, show_default=True, help="Displace atoms randomly with stdev.")
@click.option("-sv", "--scale-volume", default=1.0, type=float, show_default=True, help="Scale input volume.")
@click.option("-n","--mpi-nprocs", default=2, type=int, show_default=True, help="Number of MPI processes to run ABINIT")
@click.option("-xc", default="PBE", show_default=True, type=click.Choice(["PBE", "PBEsol", "LDA"]),
              help="XC functional.")
@add_workdir_verbose_opts
def main(filepath, nn_name,
         relax_mode, fmax, steps, optimizer,
         kppa, rattle, scale_volume, mpi_nprocs, xc,
         workdir, verbose,
         ):

    import warnings
    warnings.simplefilter("ignore")

    # Get pseudos
    repo_name = {
        "PBE": "ONCVPSP-PBE-SR-PDv0.4",
        "PBEsol": "ONCVPSP-PBEsol-SR-PDv0.4",
        "LDA": "ONCVPSP-LDA-SR-PDv0.4",
    }[xc]
    print(f"Using {repo_name=}")
    pseudos = get_repo_from_name(repo_name).get_pseudos("standard")

    # Get atoms
    if os.path.exists(filepath):
        structure = Structure.from_file(filepath)
        if abs(scale_volume - 1.0) > 0.0:
            print(f"Scaling input volume by {scale_volume=}")
            structure.scale_lattice(scale_volume * structure.lattice.volume)
        atoms = get_atoms(structure)
    else:
        raise ValueError(f"Cannot init Atoms from {filepath=}")
        #atoms = bulk(filepath)

    if rattle:
        print("Displacing atoms randomly with stdev=", rattle)
        atoms.rattle(stdev=abs(rattle), seed=42)

    prof = RelaxationProfiler(atoms, pseudos, xc, kppa, relax_mode, fmax, mpi_nprocs, steps=steps,
                              verbose=verbose, optimizer=optimizer, nn_name=nn_name)
    prof.run(workdir=workdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
