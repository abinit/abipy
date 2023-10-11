#!/usr/bin/env python
"""
Script to perform several types of calculations with ASE and ML potentials.
"""
from __future__ import annotations

import sys
import os
import json
import click
import numpy as np
import abipy.ml.aseml as aseml
import abipy.tools.cli_parsers as cli

from functools import wraps
from time import time
from abipy.core.structure import Structure
from abipy.tools.printing import print_dataframe

ASE_OPTIMIZERS = aseml.ase_optimizer_cls("__all__")


def set_default(ctx, param, filepath):
    """
    to have config file for a single command:
    Based on https://stackoverflow.com/questions/46358797/python-click-supply-arguments-and-options-from-a-configuration-file
    """
    from abipy.tools.iotools import yaml_safe_load_path
    if os.path.exists(filepath):
        print("Reading CLI options from:", filepath)
        config = yaml_safe_load_path(filepath)
        print("Config options:\n", config)
        ctx.default_map = config
    return filepath


def herald(f):
    @wraps(f)
    def wrapper(*args, **kw):
        verbose = kw.get("verbose", 0)
        #print('func:%r args:[%r, %r]' % (f.__name__, args, kw))
        #raise ValueError()
        if verbose > 3:
            print('func:%r args:[%r, %r]' % (f.__name__, args, kw))
            print(f.__doc__, end=2*"\n")
            print("Command line options:")
            print(json.dumps(kw, indent=4), end="\n")

        # Set OMP_NUM_THREADS to 1 if env var is not defined.
        num_threads = cli.fix_omp_num_threads()

        t_start = time()
        exit_code = f(*args, **kw)
        t_end = time()
        print('\n%s command completed in %2.4f sec\n' % (f.__name__, t_end - t_start))
        return exit_code

    return wrapper


def add_constraint_opts(f):
    """Add CLI options to constrain atoms."""
    def mk_cbk(type):
        def callback(ctx, param, value):
            #print(f"{param=}, {value=}")
            if value is None: return None
            return [type(s) for s in value.split()]
        return callback

    f = click.option("--fix-inds", "-fi", type=str, default=None, show_default=True,
                    callback=mk_cbk(int),
                    help='Fix atoms by indices e.g. `--fix-inds "0 1"` to fix the first two atoms.')(f)
    f = click.option("--fix-symbols", "-fs", type=str, default=None, show_default=True,
                    callback=mk_cbk(str),
                    help='Fix atoms by chemical symbols e.g. `--fix-symbols "C O"`')(f)
    return f


def add_relax_opts(f):
    """Add CLI options for structural relaxations with ASE."""
    # fmax (float): total force tolerance for relaxation convergence.
    # Here fmax is a sum of force and stress forces. Defaults to 0.1.
    f = click.option("--relax-mode", "-r", default="ions", show_default=True, type=click.Choice(["no", "ions", "cell"]),
                     help='Relaxation mode.')(f)
    f = click.option("--fmax", default=0.01, type=float, show_default=True,
                     help='Stopping criterion in eV/A')(f)
    f = click.option("--pressure", default=0.0, type=float, show_default=True, help='Scalar pressure')(f)
    f = click.option("--steps", default=500, type=int, show_default=True,
                     help="Max number of steps for ASE relaxation.")(f)
    f = click.option("--optimizer", "-o", default="BFGS", show_default=True, type=click.Choice(ASE_OPTIMIZERS),
                     help="ASE optimizer class.")(f)
    return f


def add_neb_opts(f):
    """Add CLI options for NEB calculations with ASE."""
    f = click.option("--nimages", "-n", default=14, type=click.IntRange(3, None), show_default=True,
                     help='Number of NEB images including initial/final points. Must be >= 3')(f)
    f = click.option("--relax-mode", "-r", default="ions", show_default=True, type=click.Choice(["no", "ions", "cell"]),
            help="Relax initial and final structure. Use `cell` to relax ions and cell, " +
                 "`ions` to relax atomic positions only, `no` to disable relaxation")(f)
    f = click.option("--fmax", default=0.01, type=float, show_default=True, help='Stopping criterion.')(f)
    f = click.option("--pressure", default=0.0, type=float, show_default=True, help='Scalar pressure')(f)
    f = click.option("--optimizer", "-o", default="BFGS", show_default=True, type=click.Choice(ASE_OPTIMIZERS),
                     help="ASE optimizer class.")(f)
    f = click.option("--neb-method", "-m", default="aseneb", type=click.Choice(aseml.ASENEB_METHODS),
                     show_default=True, help="ASE NEB method")(f)
    f = click.option("--climb", "-c", is_flag=True, help="Use a climbing image (default is no climbing image).")(f)
    return f


def add_nprocs_opt(f):
    """Add CLI options for multiprocessing."""
    f = click.option("-np", "--nprocs", default=-1, type=int, show_default=True,
                    help='Number of processes in multiprocessing pool. -1 to let Abipy select it automatically.')(f)
    return f


def add_workdir_verbose_opts(f):
    """Add workdir and verbose options to CLI subcommand."""
    f = click.option("--workdir", "-w", default=None, type=str,
                      help="Working directory. If not specified, a temporary directory is created.")(f)
    f = click.option('-v', '--verbose', count=True, help="Verbosity level")(f)
    return f


def add_nn_name_opt(f):
    """Add CLI options to select NN potential."""
    f = click.option("--nn-name", "-nn", default="chgnet", show_default=True,
                     type=click.Choice(aseml.CalcBuilder.ALL_NN_TYPES),
                     help='ML potential to be used')(f)
    return f


def add_nn_names_opt(f):
    """Add CLI options to select multiple NN potentials."""
    f = click.option("-nns", '--nn-names', type=str, multiple=True, show_default=True,
                    help='ML potentials to use.', default=["chgnet"])(f)
    return f


def _get_nn_names(nn_names: list[str]) -> list[str]:
    """
    Additional pre-processing of nn-names option.

    --nn-names all                   --> return all NN names.
    --nn-names all-alignn-m3gnet     --> return all NN names except alignn and m3gnet
    --nn-names chgnet-               --> return all NN names except chgnet.
    """
    if "all" in nn_names:
        # Return all possibilities.
        return aseml.CalcBuilder.ALL_NN_TYPES

    if any(n.startswith("all-") for n in nn_names):
        # --nn-names all-alignn-m3gnet     --> return all NN names except alignn and m3gnet
        assert len(nn_names) == 1
        skip_names = nn_names[0].replace("all-", "").split("-")
        return [s for s in aseml.CalcBuilder.ALL_NN_TYPES if s not in skip_names]

    if any(n.endswith("-") for n in nn_names):
        # --nn-names chgnet-               --> return all NN names except chgnet.
        skip_names = [n[:-2] for n in nn_names if n.endswith("-")]
        return [s for s in aseml.CalcBuilder.ALL_NN_TYPES if s not in skip_names]

    return nn_names


@click.group()
@click.pass_context
@click.option("--seaborn", "-sns", default=None, show_default=True,
              help='Use seaborn settings. Accept value defining context in ("paper", "notebook", "talk", "poster").')
def main(ctx, seaborn):
    """Script to perform calculations with ML potentials."""
    ctx.ensure_object(dict)

    if seaborn:
        # Activate seaborn settings for plots
        import seaborn as sns
        sns.set(context=seaborn, style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)


@main.command()
@herald
@click.pass_context
@click.argument("filepath", type=str)
@add_nn_name_opt
@add_relax_opts
@add_constraint_opts
@add_workdir_verbose_opts
@click.option('--config', default='abiml_relax.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def relax(ctx, filepath, nn_name,
          relax_mode, fmax, pressure, steps, optimizer,
          fix_inds, fix_symbols,
          workdir, verbose):
    """
    Structural relaxation with ASE and ML potential.

    Usage example:

    \b
        abiml.py.py relax FILE --fmax 0.01 -r cell --optimizer FIRE -w OUT_DIR
        abiml.py.py relax FILE --fix-inds "0 3" --fix-symbols "Si O"

    where `FILE` is any file supported by abipy/pymatgen e.g. netcdf files, Abinit input, POSCAR, xsf, etc.
    or a string such as __mp-134 to fetch the structure from the MP database.

    To change the ML potential, use e.g.:

        abiml.py.py relax -nn m3gnet [...]
    """
    if filepath.startswith("__mp-"):
        print(f"Fetching structure for mp-id {filepath[2:]} from the materials project database.")
        atoms = aseml.get_atoms(Structure.from_mpid(filepath[2:]))
    else:
        atoms = aseml.get_atoms(filepath)

    aseml.fix_atoms(atoms, fix_inds=fix_inds, fix_symbols=fix_symbols)

    ml_relaxer = aseml.MlRelaxer(atoms, relax_mode, fmax, pressure, steps, optimizer,
                                 nn_name, verbose, workdir, prefix="_abiml_relax_")

    print(ml_relaxer.to_string(verbose=verbose))
    ml_relaxer.run()
    return 0


@main.command()
@herald
@click.pass_context
@click.argument("filepath", type=str)
@add_workdir_verbose_opts
@click.option('--config', default='abiml_abinit_relax.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def abinit_relax(ctx, filepath,
                 workdir, verbose):
    """
    Interact with ABINIT in hybrid relaxation mode.
    """
    ml_relaxer = aseml.MlRelaxer.from_abinit_yaml_file(filepath)
    print(ml_relaxer.to_string(verbose=verbose))
    ml_relaxer.run()
    return 0


@main.command()
@herald
@click.pass_context
@click.argument("filepath", type=str)
@add_nn_name_opt
@click.option('--temperature', "-t", default=600, type=float, show_default=True, help='Temperature in Kelvin')
@click.option('--timestep', "-ts", default=1, type=float, show_default=True, help='Timestep in fs.')
@click.option('--steps', "-s", default=1000, type=int, show_default=True, help='Number of timesteps.')
@click.option('--loginterval', "-l", default=100, type=int, show_default=True, help='Interval for record the log.')
@click.option('--ensemble', "-e", default="nvt", show_default=True,
              type=click.Choice(["nvt", "npt", "npt_berendsen"]), help='Ensemble e.g. nvt, npt.')
@add_constraint_opts
@add_workdir_verbose_opts
@click.option('--config', default='abiml_md.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def md(ctx, filepath, nn_name,
       temperature, timestep, steps, loginterval, ensemble,
       fix_inds, fix_symbols,
       workdir, verbose):
    """
    MD simulation with ASE and ML potential.

    Usage example:

    \b
        abiml.py.py md FILE --temperature 1200 --timestep 2 --steps 5000 --workdir OUT_DIR
        abiml.py.py md FILE --fix-inds "0 3" --fix-symbols "Si O"

    where `FILE` is any file supported by abipy/pymatgen e.g. netcdf files, Abinit input, POSCAR, xsf, etc.

    To change the ML potential, use e.g.:

        abiml.py.py md -nn m3gnet [...]

    To restart a MD run, use --workdir to specify a pre-existent directory.
    """
    # See https://github.com/materialsvirtuallab/m3gnet#molecular-dynamics
    atoms = aseml.get_atoms(filepath)
    aseml.fix_atoms(atoms, fix_inds=fix_inds, fix_symbols=fix_symbols)

    ml_md = aseml.MlMd(atoms, temperature, timestep, steps, loginterval, ensemble, nn_name, verbose,
                       workdir, prefix="_abiml_md_")
    print(ml_md.to_string(verbose=verbose))
    ml_md.run()
    return 0


@main.command()
@herald
@click.pass_context
@click.argument("filepaths", nargs=2, type=str)
@add_nn_name_opt
@add_neb_opts
@add_constraint_opts
@add_workdir_verbose_opts
@click.option('--config', default='abiml_neb.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def neb(ctx, filepaths, nn_name,
        nimages, relax_mode, fmax, pressure, optimizer, neb_method, climb,
        fix_inds, fix_symbols,
        workdir, verbose
    ):
    """
    NEB calculation with ASE and ML potential.

    Usage example:

    \b
        abiml.py.py neb START_FILE END_FILE --nimages 6 --fmax=0.05 --optimizer FIRE -w OUT_DIR
        abiml.py.py neb START_FILE END_FILE --neb-method improvedtangent --climb
        abiml.py.py neb START_FILE END_FILE --fix-inds "0 3" --fix-symbols "Si O"

    where `FILE` is any file supported by abipy/pymatgen e.g. netcdf files, Abinit input, POSCAR, xsf, etc.

    To change the ML potential, use e.g.:

        abiml.py.py neb -nn m3gnet [...]
    """
    initial_atoms = aseml.get_atoms(filepaths[0])
    aseml.fix_atoms(initial_atoms, fix_inds=fix_inds, fix_symbols=fix_symbols)
    final_atoms = aseml.get_atoms(filepaths[1])
    aseml.fix_atoms(final_atoms, fix_inds=fix_inds, fix_symbols=fix_symbols)

    ml_neb = aseml.MlNeb(initial_atoms, final_atoms,
                         nimages, neb_method, climb, optimizer, relax_mode, fmax, pressure,
                         nn_name, verbose, workdir, prefix="_abiml_neb_")
    print(ml_neb.to_string(verbose=verbose))
    ml_neb.run()
    return 0


@main.command()
@herald
@click.pass_context
@click.argument("filepaths", nargs=-1, type=str) # , help="Files with structures")
@add_nn_name_opt
@add_neb_opts
@add_constraint_opts
@add_workdir_verbose_opts
@click.option('--config', default='abiml_mneb.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def mneb(ctx, filepaths, nn_name,
         nimages, relax_mode, fmax, pressure, optimizer, neb_method, climb,
         fix_inds, fix_symbols,
         workdir, verbose):
    """
    Multi-NEB calculation with ASE and ML potential.

    Usage example:

    \b
        abiml.py.py mneb FILE1 FILE2 FILE2 ... --nimages 6 --fmax=0.05 -w OUT_DIR
        abiml.py.py mneb FILE1 FILE2 FILE2 ... --neb-method improvedtangent --climb
        abiml.py.py mneb FILE1 FILE2 FILE2 ... --fix-inds "0 3" --fix-symbols "Si O"

    where `FILE` is any file supported by abipy/pymatgen e.g. netcdf files, Abinit input, POSCAR, xsf, etc.

    To change the ML potential, use e.g.:

        abiml.py.py mneb -nn m3gnet [...]
    """
    # Fix atoms
    atoms_list = [aseml.get_atoms(p) for p in filepaths]
    for atoms in atoms_list:
        aseml.fix_atoms(atoms, fix_inds=fix_inds, fix_symbols=fix_symbols)

    mneb = aseml.MultiMlNeb(atoms_list, nimages, neb_method, climb, optimizer, relax_mode, fmax, pressure,
                            nn_name, verbose, workdir, prefix="_abiml_mneb_")
    print(mneb.to_string(verbose=verbose))
    mneb.run()
    return 0


@main.command()
@herald
@click.pass_context
@click.argument("filepath", type=str)
@add_nn_names_opt
@click.option("--supercell", "-s", nargs=3, type=int, default=(2, 2, 2), show_default=True, help="Supercell dimensions.")
@click.option("--distance", "-d", type=float, show_default=True, default=0.01, help="Displacement distance in Ang.")
#@click.option('--asr', type=int, default=2, show_default=True, help="Restore the acoustic sum rule on the interatomic force constants.")
#@click.option('--dipdip', type=int, default=1, show_default=True, help="Treatment of dipole-dipole interaction.")
@click.option('--line-density', "-ld", default=20, type=float, show_default=True, help="Line density to generate the q-path for PH bands.")
@click.option('--qppa', "-qppa", default=None, type=float, show_default=True, help="q-points per atom to generate the q-mesh for PH DOS.")
@add_relax_opts
@add_workdir_verbose_opts
@click.option('--config', default='abiml_ph.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def ph(ctx, filepath, nn_names,
       supercell, distance,
       #asr, dipdip,
       line_density, qppa,
       relax_mode, fmax, pressure, steps, optimizer,
       workdir, verbose):
    """
    Use phonopy and ML potential to compute phonons.

    Usage example:

    \b
        abiml.py.py ph FILE --distance 0.03 --supercell 2 2 2

    where `FILE` provides the crystalline structure
    or a string such as __mp-134 to fetch the structure from the MP database.

    To specify the list of ML potential, use e.g.:

        abiml.py.py ddb -nn-names m3gnet --nn-names chgnet [...]

    To use all NN potentials supported, use:

        -nn-names all [...]
    """
    if filepath.startswith("__mp-"):
        print(f"Fetching structure for mp-id {filepath[2:]} from the materials project database.")
        structure = Structure.from_mpid(filepath[2:])
    else:
        structure = Structure.from_file(filepath)

    from abipy.ml.ml_phonopy import MlPhonopy
    supercell = np.eye(3) * np.array(supercell)

    nn_names = _get_nn_names(nn_names)
    ml_ph = MlPhonopy(structure, supercell,
                      distance,
                      # asr, dipdip,
                      line_density, qppa,
                      relax_mode, fmax, pressure, steps, optimizer, nn_names,
                      verbose, workdir, prefix="_abiml_ph_",
                      )
    print(ml_ph.to_string(verbose=verbose))
    ml_ph.run()
    return 0


@main.command()
@herald
@click.pass_context
@click.argument("ddb_filepath", type=str)
@add_nn_names_opt
@click.option("--supercell", "-s", nargs=3, type=int, default=(-1, -1, -1), show_default=True, help="Supercell. If < 0, supercell is taken from DDB ngqpt.")
@click.option("--distance", "-d", type=float, show_default=True, default=0.01, help="Displacement distance in Ang.")
@click.option('--asr', type=int, default=2, show_default=True, help="Restore the acoustic sum rule on the interatomic force constants.")
@click.option('--dipdip', type=int, default=1, show_default=True, help="Treatment of dipole-dipole interaction.")
@click.option('--line-density', "-ld", default=20, type=float, show_default=True, help="Line density to generate the q-path for PH bands.")
@click.option('--qppa', "-qppa", default=None, type=float, show_default=True, help="q-points per atom to generate the q-mesh for PH DOS.")
@add_relax_opts
@add_workdir_verbose_opts
@click.option('--config', default='abiml_phddb.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def phddb(ctx, ddb_filepath, nn_names,
          supercell, distance, asr, dipdip, line_density, qppa,
          relax_mode, fmax, pressure, steps, optimizer,
          workdir, verbose):
    """
    Use phonopy and ML potential to compute phonons and compare with DDB.

    Usage example:

    \b
        abiml.py.py phddb DDB_FILE --distance 0.03 --dipdip 0 --supercell 2 2 2

    where `DDB_FILE` is an Abinit DDB file
    or a string such as __mp-134 to fetch the DDB from the MP database.

    To specify the list of ML potential, use e.g.:

        abiml.py.py phddb -nn-names m3gnet --nn-names chgnet [...]

    To use all NN potentials supported, use:

        -nn-names all [...]
    """
    if ddb_filepath.startswith("__mp-"):
        print(f"Fetching DDB for mp-id {ddb_filepath[2:]} from the materials project database.")
        from abipy.dfpt.ddb import DdbFile
        with DdbFile.from_mpid(ddb_filepath[2:]) as ddb:
            ddb_filepath = ddb.filepath

    from abipy.ml.ml_phonopy import MlPhonopyWithDDB
    if any(s <= 0 for s in supercell):
        supercell = None
    else:
        supercell = np.eye(3) * np.array(supercell)

    nn_names = _get_nn_names(nn_names)
    ml_phddb = MlPhonopyWithDDB(ddb_filepath,
                                distance, asr, dipdip, line_density, qppa,
                                relax_mode, fmax, pressure, steps, optimizer, nn_names,
                                verbose, workdir, prefix="_abiml_phddb_",
                                supercell=supercell,
                                )
    print(ml_phddb.to_string(verbose=verbose))
    ml_phddb.run()
    return 0


@main.command()
@herald
@click.pass_context
@click.argument("filepath", type=str)
@add_nn_name_opt
@click.option("--max-ns", "-m", default=100, type=int, show_default=True, help='Max number of structures')
@add_relax_opts
@add_workdir_verbose_opts
@click.option('--config', default='abiml_order.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def order(ctx, filepath, nn_name,
          max_ns, relax_mode, fmax, pressure, steps, optimizer, workdir, verbose):
    """
    Generate ordered structures from CIF with partial occupancies.

    Usage example:

    \b
        abiml.py.py order FILE --max-ns 10 --relax cell

    where `FILE` is any file supported by abipy/pymatgen e.g. netcdf files, Abinit input, POSCAR, xsf, etc.

    Based on: https://matgenb.materialsvirtuallab.org/2013/01/01/Ordering-Disordered-Structures.html
    """
    ml_orderer = aseml.MlOrderer(filepath, max_ns, optimizer, relax_mode, fmax, pressure,
                                 steps, nn_name, verbose, workdir, prefix="_abiml_order_")
    print(ml_orderer.to_string(verbose=verbose))
    ml_orderer.run()
    return 0


@main.command()
@herald
@click.pass_context
@click.argument("filepath", type=str)
@add_nn_name_opt
@click.option("-isite", "--isite", required=True,
               help='Index of atom to displace or string with the chemical element to be added to input structure.')
@click.option("--mesh", type=int, default=4, show_default=True, help='Mesh size along the smallest cell size.')
@add_relax_opts
@add_nprocs_opt
@add_workdir_verbose_opts
@click.option('--config', default='abiml_scan_relax.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def scan_relax(ctx, filepath, nn_name,
               isite, mesh,
               relax_mode, fmax, pressure, steps, optimizer,
               nprocs,
               workdir, verbose
               ):
    """
    Generate 3D mesh of (nx,ny,nz) initial positions and perform multiple relaxations
    in which all atoms are fixed except the one initially placed at the mesh point.

    Usage example:

    \b
        abiml.py.py scan-relax FILE -isite 0 --mesh 4  # Move first atom in the structure
        abiml.py.py scan-relax FILE -isite H           # Add H to the structure read from FILE.

    where `FILE` is any file supported by abipy/pymatgen e.g. netcdf files, Abinit input, POSCAR, xsf, etc.

    To change the ML potential, use e.g.:

        abiml.py.py scan-relax -nn m3gnet [...]
    """
    structure = Structure.from_file(filepath)

    from abipy.ml.relax_scanner import RelaxScanner
    scanner = RelaxScanner(structure, isite, mesh, nn_name,
                           relax_mode=relax_mode, fmax=fmax, steps=steps, verbose=verbose,
                           optimizer_name=optimizer, pressure=pressure,
                           workdir=workdir, prefix="_abiml_scan_relax_")
    print(scanner)
    scanner.run(nprocs=nprocs)
    return 0


@main.command()
@herald
@click.pass_context
@click.argument('filepaths', type=str, nargs=-1)
@add_nn_names_opt
@click.option("--traj_range", type=str, show_default=True,
              help="Trajectory range e.g. `5` to select the first 5 iterations, `1:4` to select steps 1,2,3.",
              default=None)
@click.option("-e", '--exposer', default="mpl", show_default=True, type=click.Choice(["mpl", "panel"]),
              help='Plotting backend: mpl for matplotlib, panel for web-based')
@add_nprocs_opt
@add_workdir_verbose_opts
@click.option('--config', default='abiml_compare.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def compare(ctx, filepaths,
            nn_names,
            traj_range,
            exposer,
            nprocs,
            workdir, verbose
            ):
    """
    compare ab-initio energies, forces, and stresses with ml-computed ones.

    usage example:

    \b
        abiml.py.py compare FILE --nn-names matgl --nn-names chgnet

    where `FILE` can be either a _HIST.nc or a vasprun.xml FILE.
    """
    traj_range = cli.range_from_str(traj_range)
    nn_names = _get_nn_names(nn_names)
    ml_comp = aseml.MlCompareWithAbinitio(filepaths, nn_names, traj_range, verbose, workdir, prefix="_abiml_compare_")
    print(ml_comp)
    c = ml_comp.run(nprocs=nprocs)

    with_stress = True
    from abipy.tools.plotting import Exposer
    with Exposer.as_exposer(exposer, title=" ".join(os.path.basename(p) for p in filepaths)) as e:
        e(c.plot_energies(show=False))
        e(c.plot_forces(delta_mode=True, show=False))
        e(c.plot_energies_traj(delta_mode=True, show=False))
        e(c.plot_energies_traj(delta_mode=False, show=False))
        if with_stress:
            e(c.plot_stresses(delta_mode=True, show=False))
        e(c.plot_forces_traj(delta_mode=True, show=False))
        e(c.plot_stress_traj(delta_mode=True, show=False))

    return 0


@main.command()
@herald
@click.pass_context
@click.option('-v', '--verbose', count=True, help="Verbosity level")
def show_nn(ctx, verbose):
    """
    Show the NN potentials installed in the environment.
    """
    installed, versions = aseml.get_installed_nn_names(verbose=verbose, printout=True)
    return 0 if installed else 1


@main.command()
@herald
@click.pass_context
@click.option("-nns", '--nn-names', type=str, multiple=True, show_default=True,
              help='ML potentials to install.', default=["all"])
@click.option('-U', '--update', is_flag=True, default=False, show_default=True, help="Update packages.")
@click.option('-v', '--verbose', count=True, help="Verbosity level")
def install(ctx, nn_names, update, verbose):
    """
    Install NN potentials in the environment using pip
    """
    aseml.install_nn_names(nn_names=nn_names, update=update, verbose=verbose)
    installed, versions = aseml.get_installed_nn_names(verbose=verbose, printout=True)

    return 0 if installed else 1




@main.command()
@herald
@click.pass_context
@click.argument('filepaths', type=str, nargs=-1)
@add_nn_name_opt
#@add_nprocs_opt
@add_workdir_verbose_opts
@click.option('--config', default='abiml_train.yml', type=click.Path(), callback=set_default, is_eager=True, expose_value=False)
def train(ctx, filepaths,
          nn_name,
          #nprocs,
          workdir, verbose
          ):
    """
    Train a ML potential using the trajectory stored on FILE.

    usage example:

    \b
        abiml.py.py train FILE --nn-names matgl

    where `FILE` can be either a _HIST.nc or a vasprun.xml FILE.
    """
    if nn_name == "matgl":
        from abipy.ml.matgl import MatglSystem
        s = MatglSystem(filepaths, workdir, verbose)

    elif nn_name == "chgnet":
        from abipy.ml.chgnet import ChgnetSystem
        s = ChgnetSystem(filepaths, workdir, verbose)

    else:
        raise ValueError(f"Unsupported {nn_name=}")

    s.train()
    return 0


if __name__ == "__main__":
    sys.exit(main())
