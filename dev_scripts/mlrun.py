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
from abipy.ml.aseml import get_atoms, ase_optimizer_cls, CORRALGO
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
@click.option("--nn-name", "-nn", default="chgnet", show_default=True,
              type=click.Choice(["m3gnet", "matgl", "chgnet"]), help='ML potential')
@click.option("-c", "--corr-algo_str", default="delta", show_default=True,
              type=click.Choice(CORRALGO._member_names_), help='Correction algorithm')
@click.option("-algo","--algorithm", default="old", type=str, show_default=True, help="String used to test algorithms... ")

@add_relax_opts
@click.option("-k", "--kppa", default=1000, type=float, show_default=True,
              help="Defines the sampling of BZ mesh (k-points per atom)")
@click.option("--rattle", default=0.0, type=float, show_default=True, help="Displace atoms randomly with stdev.")
@click.option("-sv", "--scale-volume", default=1.0, type=float, show_default=True, help="Scale input volume.")
@click.option("-n","--mpi-nprocs", default=2, type=int, show_default=True, help="Number of MPI processes to run ABINIT")
@click.option("-xc", default="PBE", show_default=True, type=click.Choice(["PBE", "PBEsol", "LDA"]),
              help="XC functional.")
@click.option("-m","--mpi-runner", default="mpirun", type=str, show_default=True, help="String used to invoke the MPI runner. ")
@add_workdir_verbose_opts
def main(filepath, nn_name, corr_algo_str,algorithm,
         relax_mode, fmax, steps, optimizer,
         kppa, rattle, scale_volume, mpi_nprocs, xc, mpi_runner,
         workdir, verbose
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
            #print("before structure:\n", structure)
            structure = structure.scale_lattice(scale_volume * structure.lattice.volume)
            #print("after structure:\n", structure)
        atoms = get_atoms(structure)
    else:
        raise ValueError(f"Cannot init Atoms from {filepath=}")
        #atoms = bulk(filepath)

    if rattle:
        print("Displacing atoms randomly with stdev=", rattle)
        atoms.rattle(stdev=abs(rattle), seed=42)

    print("Using corr_algo:", corr_algo_str)
    corr_algo = CORRALGO.from_string(corr_algo_str)
    prof = RelaxationProfiler(atoms, pseudos, corr_algo,algorithm, xc, kppa, relax_mode, fmax, mpi_nprocs,
                              steps=steps, verbose=verbose, optimizer=optimizer, nn_name=nn_name, mpi_runner=mpi_runner)
    prof.run(workdir=workdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
